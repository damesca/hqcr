# Security audit — `hqcr` v0.2.0

**Scope.** Constant-time (CT) and zeroize audit of the pure-Rust HQC KEM
(`hqcr`), spec *HQC specifications 22/08/2025*, across both build profiles
(portable `clmul64` fallback and `+pclmulqdq`). This is a
**verification-and-documentation** pass: the KAT-verified `src/` logic is
unchanged except two small zeroize fixes (G1/G2). Two deferred limitations are
recorded as **known limitations**, not remediated.

**Method.**
- **Constant-time** — three layers: (1) manual source review tracing the
  transitive-secret chain, (2) dudect-style Welch t-test timing
  (`tests/ct_timing.rs`, `ct-audit` feature), (3) instruction-level asm + miri
  spot-checks. Each hunts three **leak classes** on secret-derived data:
  **C1** secret-dependent branch / early-exit, **C2** secret-dependent memory
  access (e.g. table index `T[secret]`), **C3** secret-dependent loop bound.
  *Timing is evidence, not proof:* a clean `|t|` does not certify CT, but a
  strong signal is positive proof of non-CT.
- **Zeroize** — per-secret lifetime traces (create → copy → drop) plus a
  compile-time `ZeroizeOnDrop` guard (`#[cfg(test)] mod zeroize_guards` in
  `lib.rs`). It checks three **defeaters**: **D1** `Copy` types (bitwise copies
  escape `ZeroizeOnDrop`), **D2** growing `Vec`/`String` (realloc leaves stale
  plaintext), **D3** `mem::forget` / leaks (skip `Drop`). Documented gaps below
  are numbered **G1–G4**.

**Threat-model boundary.** Keygen runs once on fresh randomness and is *accepted
as non-CT* (matches the HQC reference). Decapsulation runs repeatedly under a
fixed key on attacker-chosen ciphertexts — this is where CT matters.

---

## Constant-time verdict

| Location | Verdict | Notes |
|---|:--:|---|
| `mul_dense_ct` (decrypt `u·y`, secret `y`) | ✅ CT | Branch-free Karatsuba + masking `clmul64`; timing canary `\|t\|=0.90`; asm confirms mask/`cmov`, no secret `jcc` |
| `clmul64` (portable + `pclmulqdq` leaf) | ✅ CT | Portable masks per bit; SIMD = one fixed-latency `pclmulqdq`; miri-clean (crate's only `unsafe`) |
| `ct_select_key` (decaps key select) | ✅ CT | `u8::conditional_select` over all 32 B; asm = `pand`/`pandn`/`por`, no branch |
| `decaps` select / compare | ✅ CT | Decode-ok folded into `Choice`; `ct_eq` compare; only branch is the **public** ciphertext-length check |
| `rm_decode` argmax | ✅ CT | Full 128-scan, no early exit; `conditional_select`; asm = `cmovs`/`setg` + mask |
| `sample_fixed_weight_mod` (r1, r2, e) | ✅ CT | 2025-spec Barrett sampler; `ct_eq` dedup, fixed-length bit-set |
| `sample_fixed_weight` (x, y) | ⚠️ non-CT | Rejection loop branches on drawn values — **keygen-only, accepted** (reference-consistent, not attacker-queryable) |
| `mul_sparse_dense` with secret `y` | ⚠️ non-CT | Iterates set bits of `y` — **keygen-only** (`s = x + h·y`); decaps uses CT `mul_dense_ct` |
| **RS decoder + GF(2⁸) tables** | ❌ **NOT CT** | **Headline limitation** — see below |

### Headline limitation — RS decoder + GF(2⁸) tables are not constant-time

The shortened Reed-Solomon decoder (`syndromes`, `berlekamp_massey`,
`chien_search`, `forney`, `rs_decode`) and the `gf.rs` log/antilog arithmetic it
calls process the **secret-derived** decoded codeword. *"Secret" is transitive:*
in `decaps`, `m'` is secret ⇒ the decoded codeword, the RS symbols, syndromes,
locator, and **every `GF_LOG`/`GF_EXP` table index** are secret. The decoder
leaks through all three classes:

- **C1 branches on secret data** — `if s[j] != 0` / no-error early return, BM
  discrepancy and register-length branches, Chien conditional `push`, zero-operand
  branches in `gf_mul`/`gf_div`.
- **C2 secret-indexed loads** — `GF_LOG[secret]` / `GF_EXP[secret]` cache lines
  depend on the operand; `Vec` growth keyed to the secret error count.
- **C3 secret-bounded loops** — BM inner sum over the secret register length;
  Forney over a variable root count.

**Severity.** During decapsulation under a fixed key, chosen/malformed
ciphertexts induce decoding errors whose timing depends on `y` — a
**decryption-side timing oracle on the secret key**, the same class as the 2022
HQC/BIKE attacks. The implicit-rejection FO transform masks decode failure in the
*returned key value* but **not in timing**. Layer 2 measured this directly: the
isolated decoder (`codes::decode`) reads `|t| ≈ 622` (a ~9.7 % clean-vs-errored
gap) against the CT canary's `0.90` — a ~690× signal; diluted to `|t| ≈ 5.5`
through full `decaps` because the CT multiply and re-encryption dominate the cycle
count.

**Why deferred.** A CT RS/GF layer is a substantial rewrite (branch-free/bitsliced
GF multiply, fixed-iteration masked Berlekamp-Massey, no-push Chien scan,
fixed-size Forney) that risks the KAT-verified decode paths. Per the Step 19
scope (verify-and-document), it is recorded as a known limitation. Until
remediated, **`hqcr` is timing-safe for keygen/encaps but NOT hardened against a
decryption-side timing oracle.**

---

## Zeroize verdict

| Secret | Verdict | Wrapper |
|---|:--:|---|
| Secret key `(x, y)`, ephemerals `(r1, r2, e)`, ring intermediates | ✅ WIPED | `Poly<P>: #[derive(Zeroize, ZeroizeOnDrop)]` |
| All seeds (`seed_kem`, `seed_pke`, `seed_dk`), `σ`, `θ`, SHA3-512 digests | ✅ WIPED | `Zeroizing` / `DecryptionKey: ZeroizeOnDrop` |
| `m`, `m'` (encaps / recovered message) | ✅ WIPED | `Zeroizing<Vec<u8>>`, fixed length (no realloc) |
| `from_bytes` local seed | ✅ WIPED | **G1 fixed** — now wrapped in `Zeroizing` |
| decaps `k_prime` / `k_bar` candidates | ✅ WIPED | **G2 fixed** — `.zeroize()`d before return |
| RS/RM codec `Vec`s (`rs_cw`, `buf`, `poly_to_bytes`, decoder internals) | ❌ **GAP G3** | none — secret-derived, some grow (D2) |
| `expanded_secret_key_bytes` | ❌ GAP G4 | none — `#[cfg(feature="kat")]` only, negligible |
| returned `K` (`SharedKey`) | N/A | caller-owned by API design (NIST/RustCrypto convention) |
| `salt` | N/A | public (transmitted in ciphertext) |

Crate-wide defeater scan is clean: no `mem::forget`/`ManuallyDrop`/`into_raw`
(D3); no secret type is `Copy` (D1, locked by the compile-time guard); secret
`Vec`s are fixed-length (D2) — **except** the RS-codec `Vec`s.

**G3** is the same surface as the constant-time headline finding: the RS decoder
handles the secret codeword in a non-hardened way, and several `Vec`s (`sigma`
via clone, Chien `roots` via push) grow, so reallocation can leave stale
plaintext on the heap. Best hardened together with the deferred CT-decoder
rewrite to avoid two passes over KAT-verified code.

---

## Bottom line

Every security-critical path on the **encapsulation** side — ephemeral sampling,
decrypt multiply, KEM selection/compare, RM decode — is **constant-time**:
verified by source review, corroborated by timing, and confirmed at the
instruction level (`subtle` selects compile to mask/`cmov`, never a secret `jcc`;
the SIMD leaf is one `pclmulqdq`; the sole `unsafe` block is miri-clean). Every
long-term and ephemeral secret held in the KEM's own types is **zeroized on
drop**, with G1/G2 fixed.

The single material gap is the **non-constant-time Reed-Solomon decoder and its
GF(2⁸) table lookups** (and the matching un-zeroized codec `Vec`s, G3) — a
decryption-side timing oracle measured directly in Layer 2, documented as a known
limitation deferred to a future CT-decoder rewrite. No library behaviour changed
beyond the G1/G2 zeroize fixes; `cargo test` and the KAT suite pass unchanged,
and all audit instruments are gated behind the `ct-audit` feature.

> **Full detail archived.** Per-layer analysis, complete verdict tables, timing-run
> logs, and raw asm/miri dumps live on the `audit` branch under `docs/audit/`
> (`constant-time.md`, `zeroize.md`, `asm+miri_results.txt`).
