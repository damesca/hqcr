# Zeroize audit — `hqcr`

Step 19b of the roadmap (CLAUDE.md). A **verification-and-documentation** pass:
trace each secret's full lifetime and confirm a wipe on every exit path; no
KAT-verified `src/` logic is changed, and gaps are recorded as known limitations.

| | |
|---|---|
| Crate revision audited | branch `audit` |
| Method | Layer 1 manual lifetime traces + Layer 2 compile-time guard |
| Companion | `docs/audit/constant-time.md` (19a); the RS-decoder gap (G3) overlaps its §3.3 |

---

## 1. Methodology

A zeroize audit traces each secret value *create → every move / clone / copy →
drop* and confirms its backing memory is wiped on every exit path (including
panics, via `Drop`). It hunts the three classic defeaters:

- **D1 — `Copy` types.** A `Copy` value is duplicated bit-for-bit on every move;
  those copies escape any `ZeroizeOnDrop`. (Rust forbids `Copy` + `Drop`, so a
  type with a `Drop`/`ZeroizeOnDrop` impl provably cannot be `Copy`.)
- **D2 — growing `Vec`/`String`.** On reallocation the old buffer is freed
  *without* being zeroized, leaving a stale plaintext copy on the heap. Only
  fixed-size secret buffers (or wipe-before-grow) are safe.
- **D3 — `mem::forget` / leaks.** `mem::forget`, `ManuallyDrop`, `Box::leak`,
  `into_raw`, etc. skip `Drop` entirely, so no wipe runs.

Two layers:

- **Layer 1 — manual lifetime traces** *(this document)*: per-secret
  create→copy→drop with verdict + `file:line`.
- **Layer 2 — compile-time guard** *(`#[cfg(test)] mod zeroize_guards` in
  `lib.rs`)*: trait-bound assertions that fail the build if a secret type loses
  `ZeroizeOnDrop` (which also locks out `Copy`). The zeroize analog of the CT
  canary — a regression tripwire, not a runtime check.

### Crate-wide defeater scan (grep, all of `src/`)

| Defeater | Finding |
|---|---|
| **D3** `mem::forget` / `ManuallyDrop` / `Box::leak` / `into_raw` / `MaybeUninit` / `mem::*` | **No matches.** Nothing skips `Drop`. ✅ |
| **D1** `Copy` | Only `params.rs:65,69,73` — the zero-sized markers `Hqc128/192/256`, which hold no data. No secret type is `Copy`. ✅ |
| **D2** growing secret `Vec` | Secret `Vec`s (`m`, `m'`) are created at fixed length and never `push`/`extend`ed (see traces). The *only* growing secret-derived `Vec`s live in the RS decoder, which is not `Zeroizing`-wrapped at all → folded into gap **G3**. |

---

## 2. Secret inventory

Verdict legend: **WIPED** = zeroized on drop; **GAP** = secret-derived but not
wiped (documented low-risk); **N/A** = not secret.

| # | Secret | Lives as | Wrapper / wiper | Verdict | Evidence |
|--:|---|---|---|:--:|---|
| 1 | `x`, `y` (secret key) | `Poly<P>` | `#[derive(Zeroize, ZeroizeOnDrop)]` | WIPED | `poly/mod.rs:25`; created `pke.rs:109-110,169` |
| 2 | `r1`, `r2`, `e` (ephemeral) | `Poly<P>` | same | WIPED | `pke.rs:143-145` |
| 3 | ring intermediates (`hy`, `s`-pre, `uy`, `tmp`, `cm`, `sr2`, `u`, `v`) | `Poly<P>` | same (drop glue) | WIPED | `pke.rs:117-118,148-155,172-173` |
| 4 | `seed_kem` | `Zeroizing<[u8;32]>` | `Zeroizing` | WIPED | `kem.rs:65,134,154` |
| 5 | `sigma` (σ) | `Zeroizing<[u8;32]>` | `Zeroizing` (field of `DecapsulationKey`) | WIPED | `kem.rs:66,148` |
| 6 | `seed_pke` | `Zeroizing<[u8;32]>` | `Zeroizing` | WIPED | `kem.rs:146` |
| 7 | `seed_dk` | `Zeroizing<[u8;32]>` then `DecryptionKey.seed_dk` | `Zeroizing` + `DecryptionKey: ZeroizeOnDrop` | WIPED | `hash.rs:143`; `pke.rs:78-81,122` |
| 8 | `theta` (θ) | `Zeroizing<[u8;32]>` | `Theta` alias | WIPED | `hash.rs:48,104`; `pke.rs` consumes `&[u8]` |
| 9 | SHA3-512 digests in `G`/`I` (carry `K‖θ`, `seed_dk‖seed_ek`) | `Zeroizing<[u8;64]>` | `Zeroizing` | WIPED | `hash.rs:71,78` |
| 10 | `m` (encaps message) | `Zeroizing<Vec<u8>>` | `Zeroizing`, fixed length | WIPED | `kem.rs:170` |
| 11 | `m'` (decaps recovered plaintext) | `Zeroizing<Vec<u8>>` | `Zeroizing`, fixed length | WIPED | `kem.rs:221-222` |
| 12 | `DecapsulationKey` (whole) | struct | field-wise: 4,5,7 wipe; `ek` is public | WIPED | `kem.rs:64-69` |
| 13 | `from_bytes` local `seed` | `Zeroizing<[u8;32]>` | `Zeroizing` (**G1 fixed**) | WIPED | `kem.rs:88` |
| 14 | decaps `k_prime`, `k_bar` candidates | `[u8;32]` | `.zeroize()` before return (**G2 fixed**) | WIPED | `kem.rs:233-237` |
| 14b | returned shared key `K` | `[u8;32]` (`SharedKey`) | caller-owned by API design | N/A | `hash.rs:51` |
| 15 | RS/RM buffers: `rs_cw`, `buf`, `poly_to_bytes` out, decoder `Vec`s | `Vec<u8>` (plain, some grow) | — none — | **GAP G3** | `codes/mod.rs:50,54,73,99`; `codes/reed_solomon.rs` |
| 16 | `expanded_secret_key_bytes` output | `Vec<u8>` (plain) | — none — | **GAP G4** (KAT-only) | `kem.rs:116` (`#[cfg(feature="kat")]`) |
| 17 | `salt` | `[u8;16]` | — none — | N/A (public) | `kem.rs:172` |

---

## 3. Lifetime traces

### 3.1 Secret key `(x, y)` — WIPED

`keygen` (`pke.rs:99`) derives `seed_dk` (`Zeroizing`) via `i_pke_seed`, opens an
XOF on it, and samples `y` then `x` as `Poly<P>` (`pke.rs:109-110`). Both are
`ZeroizeOnDrop`; they live only to the end of `keygen` and are wiped when the
stack frame unwinds. The public `s = x + h·y` copies no secret bit out (the sum
is published). `decrypt` (`pke.rs:169`) re-derives `y` the same way, wiped at
return. No clone of `x`/`y` escapes; no `Copy`. ✅

### 3.2 Ephemerals `(r1, r2, e)` — WIPED

`encrypt` samples all three as `Poly<P>` from `XOF(θ)` (`pke.rs:143-145`); each is
`ZeroizeOnDrop`, consumed into the public `u`, `v`, and dropped (wiped) at
function exit. The intermediates `hr2`, `cm`, `sr2`, `u`, `v` are also `Poly<P>`
→ wiped. ✅

### 3.3 Seeds and σ — WIPED

`keygen_from_seed` (`kem.rs:141`) wraps the incoming `seed_kem` in `Zeroizing`
(`*seed_kem` copy, `:154`) and squeezes `seed_pke` + `sigma` as `Zeroizing`
(`:146,148`). `seed_pke` is consumed by `pke::keygen` and dropped (wiped) at the
end of `keygen_from_seed`. Inside `pke::keygen`, `i_pke_seed` returns `seed_dk`
as `Zeroizing` (`hash.rs:143`); it is copied once into `DecryptionKey.seed_dk`
(`pke.rs:122`), whose type is `ZeroizeOnDrop` — so **both** the `Zeroizing`
original and the struct field are wiped. `sigma` moves into the
`DecapsulationKey` (`Zeroizing` field, `:155→66`). The 64-byte SHA3-512 outputs
that briefly hold `K‖θ` / `seed_dk‖seed_ek` are `Zeroizing<[u8;64]>`
(`hash.rs:71`). ✅

### 3.4 θ and recovered message — WIPED

`g` returns `theta` as `Theta = Zeroizing<[u8;32]>` (`hash.rs:48,104`); the
caller passes `&theta[..]` into `encrypt` and drops (wipes) it. In `decaps`,
`m'` is `Zeroizing::new(m_prime.unwrap_or_else(|| vec![0u8; K]))`
(`kem.rs:221-222`) — fixed length, never grown (no D2), wiped on drop. ✅

### 3.5 `DecapsulationKey` — WIPED (field-wise)

`DecapsulationKey` does **not** itself derive `ZeroizeOnDrop`, and does not need
to: it owns no bare secret scalar. Its secret fields each wipe themselves —
`seed_kem`/`sigma` are `Zeroizing`, `dk_pke` is `DecryptionKey: ZeroizeOnDrop` —
and `ek` is public. Drop glue runs all three field destructors on every exit
path. `to_bytes` returns a `Zeroizing<[u8;32]>` clone (`kem.rs:78`), so even the
serialized copy is wiped. ✅

---

## 4. Documented gaps (low-risk)

**G1 — `from_bytes` local `seed` (`kem.rs:88`). ✅ FIXED.** Reconstructing a
`DecapsulationKey` copied the 32-byte `seed_KEM` (a secret — it *is* the
compressed secret key) into a plain `[u8; SEED_BYTES]` stack array that was not
wiped on return. Now wrapped in `Zeroizing` (`keygen_from_seed` only borrows it),
so it is zeroized on drop. No API or KAT change.

**G2 — decaps key candidates (`kem.rs:233-237`). ✅ FIXED (candidates).** The
intermediates `k_prime` and `k_bar` were plain `[u8; 32]` left on the stack after
the constant-time select. They are now `.zeroize()`d immediately before
`decaps` returns. The select itself is unchanged (still CT), so this adds no
behaviour and no KAT change.

> **Residual by design — returned shared key `K`.** `SharedKey = [u8; 32]` is
> `Copy` and is *returned to the caller* unwrapped; this is unavoidable (the KEM's
> purpose is to hand `K` to the caller, who then owns its lifetime — matching the
> RustCrypto / NIST KEM convention of returning a bare shared secret). It cannot
> be wrapped without a public-API change, so it stays as-is. Lifetime management
> of `K` past the call boundary is the caller's responsibility.

**G3 — RS/RM codec buffers (`codes/mod.rs`, `codes/reed_solomon.rs`).** The
concatenated codec builds and decodes through plain `Vec<u8>`: `rs_cw`/`buf` in
`encode`/`decode`, the `poly_to_bytes` output, and the RS decoder's internal
vectors (`received.to_vec()`, `corrected`, `syndromes`, `sigma`, Chien `roots`).
These hold **secret-derived** material — the encoded message and the decoded
codeword — and are not zeroized; several (`sigma` via `clone`, `roots` via
`push`) also **grow**, so reallocation can leave stale plaintext copies on the
heap (D2). *Risk:* low-to-moderate. This is the **same surface** as the headline
constant-time finding (constant-time.md §3.3): the RS decoder processes the
secret codeword in a non-hardened way. The remediation is the same deferred
CT-decoder rewrite — adding `Zeroizing` there is best done together with making
the decoder constant-time, to avoid two passes over KAT-verified code.

**G4 — `expanded_secret_key_bytes` (`kem.rs:116`, `#[cfg(feature="kat")]`).**
Emits `seed_dk‖σ‖seed_kem` as a plain `Vec<u8>` for byte-for-byte KAT comparison.
*Risk:* negligible — compiled only under the `kat` feature, never in the
production API; already documented in its own rustdoc as not-audited.

**N/A — `salt` (`kem.rs:172`).** Listed in the roadmap as "plain, dropped but not
wiped", but `salt` is **public** (transmitted in the ciphertext), so it is not a
secret and needs no wipe.

---

## 5. Layer 2 — compile-time guard

`#[cfg(test)] mod zeroize_guards` in `lib.rs` asserts the secret types satisfy
`ZeroizeOnDrop`:

```rust
fn assert_zeroize_on_drop<T: zeroize::ZeroizeOnDrop>() {}
assert_zeroize_on_drop::<Poly<Hqc128>>();   // + Hqc192, Hqc256
assert_zeroize_on_drop::<DecryptionKey>();
```

If a future refactor drops the `#[derive(ZeroizeOnDrop)]` on `Poly` or
`DecryptionKey`, the crate fails to compile. And because the derived
`ZeroizeOnDrop` carries a `Drop` impl — which Rust refuses to allow alongside
`Copy` — the same guard transitively guarantees these types can never silently
become `Copy` (defeater D1). The growing-`Vec` (D2) and forget/leak (D3)
defeaters are not expressible as a static bound; they are covered by the §1 grep
(D3: no matches) and the §3 traces (D2: secret `Vec`s are fixed-length).

---

## 6. Verdict summary

| Outcome | Items |
|---|---|
| **WIPED — verified** | 1–14 (secret key, ephemerals, ring intermediates, all seeds, σ, θ, `m`, `m'`, `DecapsulationKey`, `from_bytes` seed [G1], decaps candidates [G2]) |
| **GAP — documented low-risk** | 15 (G3 codec `Vec`s) |
| **GAP — KAT-only** | 16 (G4) |
| **N/A** | 14b (caller-owned `K`), 17 (`salt`, public) |

**Bottom line.** Every long-term and ephemeral secret that the production KEM
holds in its own types — secret key `(x, y)`, ephemerals `(r1, r2, e)`, all
seeds, σ, θ, and both copies of the recovered message — is zeroized on drop via
`Zeroize`/`ZeroizeOnDrop`/`Zeroizing`, with no `Copy`, growing-`Vec`, or
forget/leak defeater on those paths (D1/D2/D3 all clear). **G1** (the `from_bytes`
stack seed) and **G2** (the decaps `k_prime`/`k_bar` candidates) have been
**fixed** — wrapped in `Zeroizing` / `.zeroize()`d before return, with no API or
KAT change. The remaining gaps are (14b) the caller-owned returned `K`
(unavoidable without an API change) and (**G3**) the RS-codec `Vec`s — the same
secret-handling surface as the deferred constant-time RS-decoder limitation
(constant-time.md §3.3), best hardened together with that rewrite. The Layer 2
guard locks the `ZeroizeOnDrop` (and hence non-`Copy`) invariant against
regressions.
