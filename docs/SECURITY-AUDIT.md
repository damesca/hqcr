# Security audit — `hqcr` v0.3.0

**Scope.** Constant-time (CT) and zeroize audit of the pure-Rust HQC KEM
(`hqcr`), spec *HQC specifications 22/08/2025*, across both build profiles
(portable `clmul64` fallback and `+pclmulqdq`).

> **Update (Step 20 — decoder remediation).** The original audit was a
> verify-and-document pass that left one **headline limitation** open: the RS
> decoder + GF(2⁸) tables were not constant-time (a decryption-side timing
> oracle). **Step 20 (a–d) has remediated it** — branch-free GF arithmetic,
> masked Berlekamp-Massey, a branchless full-scan root finder replacing Chien
> search, and `Zeroizing` codec buffers. The verdict rows, headline section,
> zeroize **G3**, and bottom line below are updated accordingly and the timing
> was re-measured. The only remaining accepted non-CT code is **keygen-side**
> (`sample_fixed_weight`, `mul_sparse_dense`), matching the reference. Earlier
> zeroize fixes G1/G2 are retained.

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
| RS decoder + GF(2⁸) arithmetic | ✅ CT | **Remediated in Step 20** — branch-free GF, masked BM, branchless full-scan root finder, `Zeroizing` buffers; isolated `codes::decode` timing `\|t\|=0.82` (was ~622). See below |

### Headline limitation (RESOLVED in Step 20) — RS decoder + GF(2⁸) tables

**(Historical — fixed in Step 20; kept for the record and the timing before/after.)**
The shortened Reed-Solomon decoder (`syndromes`, `berlekamp_massey`,
`chien_search`, `forney`, `rs_decode`) and the `gf.rs` log/antilog arithmetic it
called processed the **secret-derived** decoded codeword. *"Secret" is transitive:*
in `decaps`, `m'` is secret ⇒ the decoded codeword, the RS symbols, syndromes,
locator, and **every `GF_LOG`/`GF_EXP` table index** were secret. The decoder
leaked through all three classes:

- **C1 branches on secret data** — `if s[j] != 0` / no-error early return, BM
  discrepancy and register-length branches, Chien conditional `push`, zero-operand
  branches in `gf_mul`/`gf_div`.
- **C2 secret-indexed loads** — `GF_LOG[secret]` / `GF_EXP[secret]` cache lines
  depend on the operand; `Vec` growth keyed to the secret error count.
- **C3 secret-bounded loops** — BM inner sum over the secret register length;
  Forney over a variable root count.

**Severity (pre-fix).** During decapsulation under a fixed key, chosen/malformed
ciphertexts induced decoding errors whose timing depended on `y` — a
**decryption-side timing oracle on the secret key**, the same class as the 2022
HQC/BIKE attacks. The implicit-rejection FO transform masks decode failure in the
*returned key value* but **not in timing**. Layer 2 measured this directly: the
isolated decoder (`codes::decode`) read `|t| ≈ 622` (a ~9.7 % clean-vs-errored
gap) against the CT canary's `0.90` — a ~690× signal; diluted to `|t| ≈ 5.5`
through full `decaps` because the CT multiply and re-encryption dominate the cycle
count.

**Resolution (Step 20).** The decoder was rewritten constant-time, following the
reference recipe below (with one deviation — a branchless full-scan root finder
instead of the additive FFT; see note):

- **20a — GF(2⁸) made branch-free / table-free.** `gf_mul` = carry-less multiply +
  fixed reduction; `gf_inv` = fixed addition chain (`a²⁵⁴`); `gf_div` branch-free.
  `GF_LOG`/`GF_EXP` are no longer indexed by secrets — kills **C2** and the
  `gf_mul`/`gf_div` zero-branch **C1**. Cross-checked exhaustively (all 256²
  products, all 255 inverses) against the old tables.
- **20b — masked Berlekamp-Massey** (port of the reference `compute_elp`): fixed
  `2δ` iterations, every update mask-merged under a `core::hint::black_box` barrier,
  inner loops to public bounds — kills the BM **C1/C3**. Cross-checked against the
  old BM on every correctable pattern.
- **20c — branchless full-scan root finder + inline Forney**, replacing Chien
  search: scans every position, tests `σ(α^{-i}) == 0` with `subtle::ct_eq`, applies
  a `conditional_select`-masked correction at the **public** index; validity folded
  to one bit. Removes the conditional `push` (**C1**) and the secret-length `Vec`
  (**G3**). All of `rs_decode` now runs a fixed sequence (no syndrome fast path, no
  early exit).
- **20d — `Zeroizing` codec buffers** (closes **G3**; see zeroize section).

**Verified.** `cargo test --features kat` passes byte-for-byte (all three sets);
the isolated decoder timing dropped from `|t| ≈ 622` to **`|t| = 0.82`**
(indistinguishable from the `mul_dense_ct` canary `|t| = 1.10`), and full `decaps`
from `|t| ≈ 5.5` to **`|t| = 1.29`** — all `< 5`, no leak signal. Dudect-style Welch
t-test, Hqc128, ~150k / 50k / 25k samples (`tests/ct_timing.rs`, `ct-audit` feature).

**Cost (correcting the earlier estimate).** Constant-time GF is intrinsically
slower than table lookups, so the decoder — formerly a small fraction of `decaps`
— now dominates the RS path. Criterion measures `decaps` **≈3.2×–4.0× slower**
than the old table-based decoder: hqc128 +225 % (≈2.5 ms), hqc192 +241 %
(≈7.2 ms), hqc256 +295 % (≈17.3 ms). **This supersedes the earlier "low
single-digit %" estimate, which was wrong** — it assumed the decoder stayed
table-fast. Keygen and encaps are unaffected (they never decode). Most of the cost
is recoverable via the deferred additive-FFT root finder + batched inversions
(Step 20c-opt) — a pure speed optimization; the current code is already
constant-time and KAT-correct.

> **Deviation note.** The reference uses a Gao–Mateer additive FFT for root
> finding; `hqcr` uses a constant-time **full evaluation** instead (same CT
> guarantee, far simpler to verify against the decode tests). The FFT is deferred
> as a performance optimization — Step 20c-opt in `CLAUDE.md`.
>
> **Verification depth.** The Step-20 decoder is verified by source review +
> dudect timing + the KAT/oracle test suite. Unlike the pre-Step-20 CT paths it
> has **not** yet had a Layer-3 asm/miri pass; recommended as follow-up.

### Remediation recipe — the official HQC reference is fully CT here

Re-audit of the official C reference (`c_implementations/reference-hqc-1`, files
titled *"Constant time implementation"*) confirms this leak is **not inherent to
HQC**: the reference closes all three classes with a concrete, KAT-blessed recipe.
This makes our gap a *fixable* limitation, not a fundamental one. The reference
techniques, mapped to our leak classes:

- **C2 — kill the GF table indices (`gf.c`).** No `GF_LOG`/`GF_EXP` on the hot
  path. `gf_mul` = branch-free carry-less multiply (`gf_carryless_mul`, mask-select
  over a 4-entry window, no zero special-case) + fixed-step `gf_reduce`.
  `gf_inverse` is a **fixed addition chain** (`a²·a³…a²⁵⁴`), not
  `GF_EXP[255 − GF_LOG[a]]`. The residual `gf_exp`/`gf_log` tables are indexed
  **only by public** loop counters / the fixed FFT basis, never by a secret.
- **C1 — kill the Berlekamp branches (`reed_solomon.c:compute_elp`).** Constant-time
  Berlekamp-Massey: every conditional update is mask-merged (`mask1`/`mask2`/`mask12`),
  with a `volatile` barrier to stop the compiler re-introducing a branch; loop
  bound is the public `2·δ`.
- **C1/C3 — replace Chien search with an additive FFT (`fft.c`).** Root finding
  evaluates σ at **every** field element via a Gao–Mateer / Bernstein–Chou–Schwabe
  additive FFT (`fft` + `fft_retrieve_error_poly`), reading off roots branchlessly
  (`error[index] ^= 1 ^ (−w[i] >> 15)`). No conditional `push`, no early exit, no
  secret-bounded loop. Error values (`compute_error_values`) are likewise computed
  with masked, full-length scans and a branch-free position counter.
- **Syndromes** use a precomputed `alpha_ij_pow[2δ][n1−1]` table indexed only by
  public counters.

(The RM decoder, `reed_muller.c`, is already CT in both the reference and `hqcr` —
full 128-scan FHT + branchless `find_peaks` argmax.)

Porting plan and its performance implications are written up as **Step 20** in
`CLAUDE.md`.

---

## Zeroize verdict

| Secret | Verdict | Wrapper |
|---|:--:|---|
| Secret key `(x, y)`, ephemerals `(r1, r2, e)`, ring intermediates | ✅ WIPED | `Poly<P>: #[derive(Zeroize, ZeroizeOnDrop)]` |
| All seeds (`seed_kem`, `seed_pke`, `seed_dk`), `σ`, `θ`, SHA3-512 digests | ✅ WIPED | `Zeroizing` / `DecryptionKey: ZeroizeOnDrop` |
| `m`, `m'` (encaps / recovered message) | ✅ WIPED | `Zeroizing<Vec<u8>>`, fixed length (no realloc) |
| `from_bytes` local seed | ✅ WIPED | **G1 fixed** — now wrapped in `Zeroizing` |
| decaps `k_prime` / `k_bar` candidates | ✅ WIPED | **G2 fixed** — `.zeroize()`d before return |
| RS/RM codec `Vec`s (`rs_cw`, `buf`, `poly_to_bytes`, decoder internals) | ✅ WIPED | **G3 fixed (Step 20d)** — all `Zeroizing`; fixed-size after the 20c rewrite (no realloc) |
| `expanded_secret_key_bytes` | ❌ GAP G4 | none — `#[cfg(feature="kat")]` only, negligible |
| returned `K` (`SharedKey`) | N/A | caller-owned by API design (NIST/RustCrypto convention) |
| `salt` | N/A | public (transmitted in ciphertext) |

Crate-wide defeater scan is clean: no `mem::forget`/`ManuallyDrop`/`into_raw`
(D3); no secret type is `Copy` (D1, locked by the compile-time guard); secret
`Vec`s are fixed-length (D2) — **including** the RS-codec `Vec`s after Step 20d.

**G3 (fixed in Step 20d).** Formerly the RS decoder handled the secret codeword
through un-wiped `Vec`s, some of which grew (`sigma` via `clone`, Chien `roots`
via `push`) so reallocation could leave stale plaintext on the heap. The 20c
rewrite removed the growing `Vec`s (fixed-size buffers, no `clone`/`push`), and
20d wraps every secret-derived codec buffer — `s`, `σ`, `Ω`, `σ'`, `corrected`,
`s2`, the BM scratch (`x_sigma_p`, `sigma_copy`), `rs_cw`/`buf` in `codes/mod.rs`,
and the `rs_encode` dividend `d` — in `Zeroizing`. Buffers that are *returned*
(`syndromes`→`s`, BM→`σ`, `error_evaluator`→`Ω`, `formal_derivative`→`σ'`,
`poly_to_bytes`→`buf`) are wrapped at the call site, so the same allocation is
wiped on the caller's drop. **D2 no longer applies.** (RM decoding uses fixed
stack arrays, not heap `Vec`s, so it was never a realloc surface.)

---

## Bottom line

Every security-critical path under a fixed key — ephemeral sampling, decrypt
multiply (`u·y`), the **RS/GF decoder** (Step 20), KEM selection/compare, RM
decode — is now **constant-time**: verified by source review, corroborated by
timing (isolated decoder `|t| = 0.82`, full `decaps` `|t| = 1.29`, vs the
`mul_dense_ct` canary `|t| = 1.10`; all `< 5`), and — for the pre-Step-20 paths —
confirmed at the instruction level (`subtle` selects compile to mask/`cmov`, never
a secret `jcc`; the SIMD leaf is one `pclmulqdq`; the sole `unsafe` block is
miri-clean). Every long-term and ephemeral secret, **including all RS/RM codec
buffers** (G3, Step 20d), is zeroized on drop (G1/G2 also fixed).

The previously headline **decryption-side timing oracle (RS decoder + GF tables)
is closed.** The remaining accepted non-CT code is **keygen-only** —
`sample_fixed_weight` (rejection sampling) and `mul_sparse_dense` (iterates secret
bit positions) — which runs once on fresh randomness, is not attacker-queryable,
and matches the HQC reference. **Cost:** `decaps` is ≈3.2×–4.0× slower than the old
table-based decoder (see Headline § *Cost*); a deferred additive-FFT root finder
(Step 20c-opt) would recover most of it. Open follow-ups: a Layer-3 asm/miri pass
over the Step-20 decoder, and **G4** (`expanded_secret_key_bytes`, KAT-only).
`cargo test` and the KAT suite pass unchanged; all audit instruments are gated
behind the `ct-audit` feature.

> **Full detail archived.** Per-layer analysis, complete verdict tables, timing-run
> logs, and raw asm/miri dumps live on the `audit` branch under `docs/audit/`
> (`constant-time.md`, `zeroize.md`, `asm+miri_results.txt`).
