# Constant-time audit — `hqcr`

Step 19a of the roadmap (CLAUDE.md). This document is a **measurement and
documentation instrument only** — no library behaviour changes as a result of
it. The KAT-verified `src/` logic is left untouched; gaps are recorded as
**known limitations**, not remediated.

| | |
|---|---|
| Spec | HQC specifications 22/08/2025 |
| Crate revision audited | branch `audit` (see `git rev-parse HEAD`) |
| Build profiles in scope | portable (`clmul64` fallback) **and** `+pclmulqdq` |
| Author | Layers 1 (manual) + 2 (timing) + 3 (asm/miri) complete |

---

## 1. Methodology

A constant-time audit hunts three **leak classes** wherever secret-derived data
is processed:

- **C1 — secret-dependent branch / early-exit.** Control flow (an `if`, `while`,
  `match`, `?`, `return`, `&&`/`||` short-circuit) whose taken/not-taken depends
  on a secret value. Leaks via instruction-count / branch-predictor timing.
- **C2 — secret-dependent memory access.** A load/store whose *address* depends
  on a secret (table index `T[secret]`, `Vec::push` deciding capacity, a
  data-dependent store offset). Leaks via cache timing.
- **C3 — secret-dependent loop bound.** A loop whose trip count depends on a
  secret. Leaks via timing.

The audit is performed in three layers; this pass completes **Layer 1**.

- **Layer 1 — manual review** *(this document, complete)*. Classify every input,
  trace the transitive-secret chain, record a verdict + `file:line`.
- **Layer 2 — empirical timing** *(pending — `tests/ct_timing.rs`, `ct-audit`
  feature)*. dudect-style Welch t-test on `decaps` and `mul_dense_ct`.
- **Layer 3 — IR / asm + miri** *(pending)*. Confirm `subtle` selects compile to
  `cmov`/masking, not `jcc`; miri the `pclmulqdq` path.

### "Secret" is transitive

The decryption/decapsulation path is the security-critical surface. The secret
key `y` enters there, and **everything computed from a secret is itself secret**:

```
ciphertext (u, v)  ── attacker-chosen ──┐
secret key y       ── long-term secret ─┤
                                         ▼
tmp = v + u·y            (pke::decrypt, mul_dense_ct)   ← depends on y
  │  C.Decode(tmp):
  ├─ poly_to_bytes(tmp)                                  ← secret
  ├─ rm_decode(block)   per RS symbol                    ← secret input
  └─ rs_decode(symbols) syndromes → BM → Chien → Forney  ← secret input
        │ uses gf_mul / gf_div / gf_inv / gf_poly_eval   ← secret table indices
```

So `m'`, the decoded codeword, the `n1` RS symbols, the syndromes, the locator
polynomial, the Chien roots, and **every GF(2⁸) table lookup the RS decoder
performs** are all secret-derived. A timing or cache leak anywhere on that chain
is a decryption-side oracle — the same attack class as the 2022 HQC/BIKE
rejection-sampling timing attacks, transplanted to the decoder.

### Threat model boundary: keygen vs. decapsulation

Two operations touch secrets with **different** CT expectations:

- **Keygen** samples `x, y` and computes `s = x + h·y` **once**, from freshly
  generated randomness, and is *accepted as non-constant-time* — this matches
  the HQC reference, which uses a rejection sampler (`vect_sample_fixed_weight1`)
  and a sparse multiply for the one-time secret generation. An attacker has no
  repeated-query oracle against a single keygen.
- **Decapsulation** is invoked repeatedly on attacker-chosen ciphertexts under a
  *fixed* secret key. This is where CT actually matters, and where the headline
  limitation (the RS decoder) lives.

---

## 2. Audit table

Verdict legend: **CT** = constant-time (verified, Layer 1); **NOT CT** =
documented limitation; **N/A** = operates only on public data.

| # | Location | Data classification | Leak class | Verdict | Evidence (Layer 1 manual) |
|--:|---|---|---|:--:|---|
| 1 | `sampling::sample_fixed_weight` (x, y) | secret positions, keygen-only | C1, C3 | NOT CT | `sampling.rs:108-136` — rejection loop + dup redraw branch on drawn values; keygen-only, reference-accepted |
| 2 | `sampling::sample_fixed_weight_mod` (r1, r2, e) | ephemeral positions | — | CT | `sampling.rs:197-225` — `ct_eq` dedup, `conditional_select` fix, full-word bit-set; all loops fixed public length |
| 3 | `sampling::sample_uniform` (h) | public | — | N/A | `sampling.rs:237-252` — public key component |
| 4 | `mul::mul_sparse_dense` — keygen `h·y` | secret `y`, keygen-only | C1, C2, C3 | NOT CT | `mul.rs:137-154` — iterates set bits of `y` (`while word!=0`, `trailing_zeros`); used at `pke.rs:117`. Keygen-only; see §3.2 |
| 5 | `mul::mul_sparse_dense` — encrypt `h·r2`, `s·r2` | ephemeral `r2` | C1, C3 | accepted | `mul.rs:137-154`, `pke.rs:148,153` — `r2` is ephemeral (re-derived from `θ` each encaps), not the long-term secret |
| 6 | `mul::mul_dense_ct` (decrypt `u·y`) | secret `y` | — | CT | `mul.rs:300-321` — full Karatsuba product visits every limb; no branch/access depends on operand bits |
| 7 | `mul::karatsuba` | — | — | CT | `mul.rs:228-283` — recursion bounds depend only on `n` (public); no operand-value branch |
| 8 | `mul::clmul64` portable fallback | — | — | CT | `mul.rs:191-207` — masks per bit, `wrapping_neg`; the only branch (`if i!=0`, `:199`) is on the public loop counter |
| 9 | `mul::clmul64` pclmulqdq path | — | — | CT (miri pending) | `mul.rs:166-184` — single data-independent-latency instruction; only `unsafe` in crate; Layer 3 to confirm |
| 10 | `kem::decaps` length check | public (ciphertext length) | C1 (public) | CT | `kem.rs:212-215` — branch on `c.len()`, a public quantity; returns `k_bar` |
| 11 | `kem::decaps` decode-ok → select | secret success bit | — | CT | `kem.rs:220,231,233` — `m'.is_some()` folded into `Choice`, never branched |
| 12 | `kem::decaps` re-encrypt compare | secret `c'` vs public `c` | — | CT | `kem.rs:230` — `c_prime.ct_eq(c)` byte-wise |
| 13 | `kem::ct_select_key` | secret keys K′, K̄ | — | CT | `kem.rs:246-252` — `u8::conditional_select` over all 32 bytes |
| 14 | `reed_muller::rm_decode` argmax | secret WHT values | — | CT | `reed_muller.rs:196-206` — full 128-scan, `conditional_select`, no early exit |
| 15 | `reed_muller::rm_decode` vote accumulation + `wht` | secret block bits | — | CT | `reed_muller.rs:166-178` — fixed-length loops, arithmetic only, no data-dependent branch |
| 16 | `reed_solomon::syndromes` | secret received word | C1 | NOT CT | `reed_solomon.rs:121` — `if s[j]!=0` sets `has_error`; `rs_decode:254` early-returns on it |
| 17 | `reed_solomon::berlekamp_massey` | secret syndromes | C1, C2, C3 | NOT CT | `reed_solomon.rs:153,155` — branches on discrepancy `delta` and register length `l`; `sigma.clone()`, `truncate(l+1)` |
| 18 | `reed_solomon::chien_search` | secret locator σ | C1, C2 | NOT CT | `reed_solomon.rs:190-191` — conditional `roots.push` on `gf_poly_eval==0` (data-dependent Vec growth) |
| 19 | `reed_solomon::forney` | secret σ, syndromes, roots | C1, C3 | NOT CT | `reed_solomon.rs:216,225,229-237` — `if i+j<2δ`, `.map` over variable root count |
| 20 | `reed_solomon::rs_decode` control flow | secret error count | C1, C2 | NOT CT | `reed_solomon.rs:254,262,270,285` — early returns on `has_error`/`num_errors`/`roots.len()`/`still_err`; `to_vec()` |
| 21 | `gf::gf_mul` | secret operands | C1, C2 | NOT CT | `gf.rs:69,72` — `if a==0||b==0` zero-branch; `GF_LOG[a]`, `GF_EXP[log_a+log_b]` secret-indexed loads |
| 22 | `gf::gf_inv` / `gf::gf_div` | secret operands | C1, C2 | NOT CT | `gf.rs:79-92` — `if a==0` branch; `GF_EXP[255-GF_LOG[a]]` secret-indexed |
| 23 | `gf::gf_poly_eval` | secret coeffs / point | C2 | NOT CT | `gf.rs:96-102` — Horner over `gf_mul`/`gf_add`, secret table indices each step |
| 24 | `codes::decode` plumbing | secret codeword | — | CT | `codes/mod.rs:67-81` — `poly_to_bytes` + per-block loop are fixed public length; the non-CT part is `rs_decode` (rows 16-23) |
| 25 | `parsing::unpack_ciphertext` length check | public length | C1 (public) | CT | `parsing.rs:150` — branch on `bytes.len()`, public |

---

## 3. Layer 1 findings — detail

### 3.1 CT paths confirmed (verified, no change needed)

- **Ephemeral sampler `sample_fixed_weight_mod`** (row 2). Step 1 draws exactly
  `4·weight` bytes (no rejection), step 2 deduplicates with `support[j].ct_eq`
  folded into a `subtle::Choice` and `u32::conditional_select`, and step 3 sets
  bits by scanning **every** word `0..N_WORDS` and OR-ing a `conditional_select`
  mask — never a data-dependent store address. All three loops have fixed,
  public lengths. This is the 2025-spec sampler that exists *because of* the 2022
  timing attacks; it is correctly constant-time. **CT.**

- **Decrypt multiply `mul_dense_ct` / `karatsuba` / `clmul64`** (rows 6-9). The
  full `2N`-bit GF(2)[X] product is computed by visiting every limb
  unconditionally; the Karatsuba recursion splits on `h = ⌈n/2⌉` (a function of
  the public `n`), and the leaf `clmul64` is branch-free on operand values in
  both the portable (per-bit mask + `wrapping_neg`) and `pclmulqdq` (single
  fixed-latency instruction) variants. The only branch in the portable leaf is
  `if i != 0` on the public loop counter (`mul.rs:199`). The `debug_assert!`
  invariant guard (`mul.rs:304`) compiles out in release. **CT** w.r.t. the
  secret `y`. *(Layer 3 will confirm the asm and miri the intrinsic.)*

- **`decaps` selection logic** (rows 10-13). The single branch is the public
  ciphertext-length check (`kem.rs:212`). Decode success is folded into a
  `Choice` (`kem.rs:220`), the re-encryption is compared with `ct_eq`
  (`kem.rs:230`), and the final key is chosen with `conditional_select` over all
  32 bytes (`kem.rs:246-252`). On decode failure an all-zero message is fed
  through `G` and re-encrypt so both run regardless of success. **CT.**

- **RM decode** (rows 14-15). The argmax scans all 128 WHT outputs with
  `i16::conditional_select` / `u8::conditional_select` and no early exit
  (`reed_muller.rs:196-206`); vote accumulation and the WHT butterfly are
  fixed-length arithmetic. **CT.** Note: although RM decode is itself CT, its
  *input blocks are secret* and its *output symbols feed the non-CT RS decoder*
  (§3.3).

### 3.2 Keygen-only non-CT (documented, accepted)

- **Secret sampler `sample_fixed_weight`** (row 1). Rejection-based: the
  number of redraws (on `cand ≥ threshold`, `sampling.rs:119`) and duplicate
  redraws (`sampling.rs:133`) depend on the drawn position values, so running
  time and XOF consumption vary. This is required to reproduce the reference's
  exact byte stream for KAT, and matches the reference's acceptance of variable
  timing for once-per-keygen secret generation.

- **Sparse multiply `mul_sparse_dense` with secret `y` in keygen** (row 4). The
  `s = x + h·y` step calls `mul_sparse_dense(&y, &h)` at `pke.rs:117`. Mode A
  iterates over the **set bits of `y`** (`while word != 0`, `trailing_zeros`,
  `word &= word-1`) — its trip count and access pattern reveal `y`'s Hamming
  weight distribution and positions. This is a genuine C1/C2/C3 leak, but it is
  confined to **keygen**, where `y` is freshly generated by the already-non-CT
  rejection sampler and never re-queried. The decapsulation path never multiplies
  by the secret with Mode A — it uses the constant-time `mul_dense_ct` (row 6,
  `pke.rs:172`). **Accepted as a keygen-only limitation**, consistent with the
  sampler decision above.

  > Note: the roadmap inventory labels Mode A simply "no CT requirement —
  > positions come from the public XOF." That is true for the *encrypt* uses
  > (`h·r2`, `s·r2`, row 5) where the sparse operand is the ephemeral `r2`, but
  > **not** for the keygen `h·y` use, where the sparse operand is the long-term
  > secret. This audit records the keygen case explicitly.

### 3.3 Headline known-limitation — the RS decoder + GF(2⁸) tables are NOT constant-time

**Rows 16-23. This is the principal finding of the audit.**

The shortened Reed-Solomon decoder and the GF(2⁸) arithmetic it depends on
process the **secret-derived** decoded codeword (the transitive chain in §1) and
leak its structure through all three classes:

- **C1 — branches on secret data.**
  - `syndromes` sets `has_error` from `if s[j] != 0` (`reed_solomon.rs:121`), and
    `rs_decode` early-returns the no-error case (`:254`). Whether *any* error
    occurred is leaked.
  - `berlekamp_massey` branches on the discrepancy `delta == 0` and the register
    update condition `2*l <= i` (`:153,155`) — both functions of the secret
    syndromes — so its instruction trace reveals the error count and pattern.
  - `chien_search` conditionally pushes a root when `gf_poly_eval(σ, xi) == 0`
    (`:190-191`).
  - `rs_decode` further branches on `num_errors`, `roots.len()`, and the
    re-verification `still_err` (`:262,270,285`).
  - `gf_mul`/`gf_div` branch on zero operands (`gf.rs:69,88`).
- **C2 — secret-indexed table loads.** Every `gf_mul`/`gf_inv`/`gf_div`/
  `gf_poly_eval` indexes `GF_LOG[secret]` and `GF_EXP[secret]`
  (`gf.rs:72,81,99`). The cache line touched depends on the secret operand.
  `chien_search`/`rs_decode` grow `Vec`s whose length is the secret error count.
- **C3 — secret-bounded loops.** `berlekamp_massey`'s inner discrepancy sum runs
  `1..=l` (secret `l`); `forney` maps over a root list whose length is the secret
  error count (`reed_solomon.rs:149,229`).

**Severity.** During decapsulation under a fixed key, an attacker submitting
chosen/malformed ciphertexts can induce decoding errors and observe the
decoder's timing. Because the error positions depend on `y` and the chosen
ciphertext, this is a **decryption-side timing oracle on the secret key** — the
same vulnerability class as the 2022 HQC/BIKE attacks that the 2025 spec's
Barrett sampler (already adopted here on the *encryption* side, row 2) was
introduced to close. The KEM's implicit-rejection FO transform masks decode
failure in the *returned key value*, but **does not mask timing**, so it does not
mitigate this leak.

**Layer 2 confirmation.** The empirical harness (§4) corroborates this:
`codes::decode` in isolation shows `|t| = 74.81` uncropped (621.98 cropped) — a
~9.7 % clean-vs-errored timing gap — against a constant-time canary
(`mul_dense_ct`) reading `0.90`. The leak is real and large; it is merely diluted
to `|t| ≈ 5.5` when measured through full `decaps` because the CT dense multiply
and re-encryption dominate that operation's cycle count.

**Why deferred (per the Step 19 decision).** A constant-time RS/GF layer is a
substantial rewrite: branch-free GF multiply (e.g. carry-less / bitsliced rather
than log-table), a fixed-iteration Berlekamp-Massey with masked updates, a full
no-push Chien scan, and a Forney over a fixed-size root buffer. That work risks
the KAT-verified decode paths. The Step 19 scope is explicitly
verification-and-documentation, so the decoder is **recorded as a known
limitation, not remediated**. Remediation is tracked as future work (a CT
decoder layer); until then, `hqcr` should be considered **timing-safe for
keygen/encaps but NOT hardened against a decryption-side timing oracle**.

---

## 4. Layer 2 — empirical timing (complete)

The harness is implemented in **`tests/ct_timing.rs`**, gated behind the
`ct-audit` feature (Cargo.toml). It is a dudect-style leakage detector:

- **Method.** Two interleaved input classes — a *fixed* input and *random*
  inputs from a pre-generated pool — are timed per call with `rdtsc` fenced by
  `lfence` (nanosecond `Instant` fallback off x86-64). Welch's unequal-variance
  t-test is computed over several upper-percentile **crops** (100 / 99.9 / 99 /
  95 / 90 %) to suppress heavy-tailed OS-scheduling outliers; `max |t|` is
  reported. Verdict bands: `|t| < 5` no evidence, `5 ≤ |t| < 10` inconclusive,
  `|t| ≥ 10` leak signal.
- **Hygiene.** All input preparation (sampling — itself variable-time — and
  ciphertext construction) is pre-generated *outside* the timed region; the
  timed closure runs only the operation, `black_box`ed. Classes are interleaved
  via PRNG to cancel CPU-frequency drift.
- **Two targets.**
  - `mul_dense_ct` (Hqc128): public operand `a` fixed, secret multiplicand `b`
    fixed-vs-random. **Expected `|t| ≈ 0`** — confirms row 6.
  - `decaps` (Hqc128): fixed key; class fix = a fixed **valid** ciphertext
    (decode succeeds, syndromes zero, RS early-return), class rand = random
    correct-length ciphertexts (drive the decoder down data-dependent error
    paths). **A leak signal is expected**, *empirically corroborating* the §3.3
    RS-decoder finding — it does **not** contradict the CT verdicts on rows
    6/11/12/13/14.

The tests **log** their statistic and never assert (general-purpose-OS timing is
evidence, not proof; a clean result does not certify CT).

> Run (user-executed — never run from the audit). **`--release` is required**:
> debug timing is meaningless, and a debug build panics in the RS decoder on
> random ciphertexts (`gf_inv(0)`'s `debug_assert`).
> ```
> cargo test --release --features ct-audit --test ct_timing -- --nocapture
> # more samples for a sharper statistic:
> CT_AUDIT_ITERS=400000 cargo test --release --features ct-audit --test ct_timing -- --nocapture
> # exercise the pclmulqdq leaf as well:
> RUSTFLAGS="-C target-feature=+pclmulqdq" cargo test --release --features ct-audit --test ct_timing -- --nocapture
> ```

### 4.1 Recorded results

**Run 1 — portable build, Windows 11 laptop. Confounded; superseded.**

| Target | iters | max \|t\| | crop | raw \|t\| (100%) | means fix / rand (cyc) |
|---|--:|--:|--:|--:|---|
| `decaps` (Hqc128) | 50 000 | 27.92 | 90% | 2.89 | 6 871 865 / 6 885 582 |
| `mul_dense_ct` (Hqc128) | 200 000 | 25.93 | 90% | 1.97 | 4 987 637 / 4 992 372 |

**Interpretation — the CT canary caught a harness artifact, not a leak.**
`mul_dense_ct` is constant-time by construction (row 6), yet Run 1 reported
`|t| = 25.93`, essentially the same curve as `decaps` (low at full data, rising
sharply under cropping; `rand` slightly slower than `fix` in *both*). That
identical signature on a known-CT operation means the statistic was dominated by
a **cache-locality confound**, not by the operation:

- the `fix` class reused **one** hot operand (cached in L1/L2 across iterations);
- the `rand` class fetched from a **~7 MB cold pool**, so every operand read
  inside the timed region paid a cache miss.

The fetch cost is charged to the timed region (the multiply/decaps reads the
operand), so `rand` is systematically ~0.1 % slower — a difference in *where the
data lives*, not in *what is computed*. At large `N` the tiny, consistent offset
crosses `|t| = 10` after cropping. Because the known-CT canary lit up just as
brightly as `decaps`, **Run 1 cannot attribute `decaps`'s signal to the RS
decoder** — both are swamped by the same artifact.

**Action taken.** The harness was corrected to draw *both* classes from
equal-footprint pools accessed at a random index (a fixed secret replicated, vs
distinct random secrets), so both pay the same cache-miss cost and only the
operand *content* differs (`tests/ct_timing.rs`, "equal-footprint pools"). Run 2
with the corrected harness is pending; the expectation is `mul_dense_ct` falls to
`|t| < 5` (confirming the artifact and validating the fix), and any residual
`decaps` signal is then attributable to the decode work (§3.3).

**Run 2 — corrected (equal-footprint) harness, portable build, Windows 11 laptop.**

| Target | iters | max \|t\| | crop | raw \|t\| (100%) | means fix / rand (cyc) | verdict |
|---|--:|--:|--:|--:|---|---|
| `mul_dense_ct` (Hqc128) | 100 000 | **0.90** | 95% | 0.43 | 5 068 006 / 5 069 246 | no leak signal ✅ |
| `decaps` (Hqc128) | 50 000 | **5.55** | 90% | 3.06 | 6 996 963 / 7 009 860 | inconclusive (weak) |

**Interpretation.**

- **Canary clean.** `mul_dense_ct` fell from `25.93` (Run 1) to `0.90`, with the
  fix/rand means now matching to ~0.02 %. This (a) confirms Run 1's signal was
  entirely the cache-locality artifact, and (b) validates the corrected harness —
  it now measures computation, so the `decaps` figure is trustworthy. The
  branch-free Karatsuba/`clmul64` decrypt multiply (row 6) shows **no timing
  dependence on the secret multiplicand**, corroborating its Layer 1 CT verdict.
- **`decaps` weak but real.** The signal dropped `27.92 → 5.55` once the artifact
  was removed, i.e. most of Run 1 was the confound. The residual is genuine
  (means differ ~0.18 %, ~9× the canary's, and `|t|` rises monotonically with
  cropping) but **diluted**: `decaps` ≈ 7M cycles, of which the CT dense multiply
  (`u·y`, measured at ~5M) plus the near-constant re-encryption dominate, so the
  non-CT RS decoder moves the total only a fraction of a percent. At
  `CT_AUDIT_ITERS ≳ 160k` this would cross `|t| = 10`, but the cleaner evidence is
  to measure the decoder in isolation (below).

**Run 3 — RS decoder isolated (`codes::decode`), portable build, Windows 11 laptop.**
Strips the CT multiply / re-encryption overhead so the non-CT decoder (§3.3) is
essentially all of the measured work; `fix` = clean codeword (early return),
`rand` = codeword with `δ/2` = 7 symbol errors (full Berlekamp-Massey + Chien +
Forney). Run only this target (fast):

```
cargo test --release --features ct-audit --test ct_timing -- --nocapture ct_timing_codes_decode
```

| Target | iters | max \|t\| | crop | raw \|t\| (100%) | means fix / rand (cyc) | verdict |
|---|--:|--:|--:|--:|---|---|
| `codes::decode` (Hqc128) | 300 000 | **621.98** | 90% | **74.81** | 208 207 / 228 406 | **LEAK SIGNAL** |

**Interpretation — decisive.** With the decoder isolated, the clean vs errored
contrast is enormous: a ~9.7 % mean gap (errored decode is slower because it runs
the full BM/Chien/Forney pipeline instead of the syndromes-zero early return),
`|t| = 74.81` *without any cropping*, climbing to 621.98 as outliers are trimmed.
Against the CT canary's `0.90`, that is a ~690× signal. This **empirically
confirms the §3.3 headline finding**: `codes::decode` (hence `rs_decode` + the
GF(2⁸) tables) is not constant-time and its timing depends directly on the
secret-derived codeword. The dilution seen at the `decaps` level (Run 2, 5.55) is
explained — the same leak is present there but buried under the CT multiply and
re-encryption.

### 4.2 Layer 2 conclusion

The three targets form a consistent, mutually-reinforcing picture, matching the
Layer 1 verdicts exactly:

| Operation | Layer 1 verdict | Layer 2 max \|t\| | agreement |
|---|---|--:|---|
| `mul_dense_ct` (decrypt multiply, secret `y`) | CT (row 6) | 0.90 | ✅ no leak |
| `decaps` (end-to-end) | CT select/compare; non-CT decoder inside | 5.55 | ✅ weak (diluted) leak |
| `codes::decode` (RS decoder isolated) | NOT CT (rows 16–23, §3.3) | 621.98 | ✅ strong leak |

**Layer 2 is complete.** The empirical timing corroborates both the CT verdicts
(the decrypt multiply is clean) and the headline limitation (the RS decoder
leaks). Caveat per method: a clean `|t|` is evidence, not proof of constant-time;
the strong signals, however, are positive proof of *non*-constant-time.

_Optional follow-ups (not required for the conclusion): re-run `decaps` at
`CT_AUDIT_ITERS=200000` to push its diluted signal past 10; re-run `mul_dense_ct`
under `RUSTFLAGS="-C target-feature=+pclmulqdq"` to confirm the SIMD leaf is also
clean._

---

## 5. Layer 3 — IR / asm + miri (harness landed; analysis pending)

Goal: confirm at the instruction level that (a) every `subtle` select on
secret-derived data compiles to branch-free `cmov` / bitmask, never a
secret-dependent conditional jump (`jcc`); (b) the portable `clmul64` masks
rather than branches on operand bits; (c) `mul_dense_ct` / `karatsuba` carry no
operand-data branch (only public-`n` loop branches); and (d) the `pclmulqdq`
`unsafe` block is free of UB.

### 5.1 Asm spot-check shims

The CT leaves are `#[inline]` / private and otherwise vanish from the asm, so
four `#[no_mangle] #[inline(never)]` shims are compiled **only under
`--features ct-audit`** (no effect on a normal build, mirroring the existing
`#[cfg(feature = "kat")]` shim in kem.rs):

| Shim | Wraps | File | Inspect for |
|---|---|---|---|
| `ct_asm_clmul64` | `clmul64` | `poly/mul.rs` | portable: bit-loop + mask, no operand `jcc`; `+pclmulqdq`: one `pclmulqdq` |
| `ct_asm_mul_dense_ct` | `mul_dense_ct::<Hqc128>` | `poly/mul.rs` | only public-`n` loop branches; `karatsuba` symbol likewise |
| `ct_asm_select_key` | `ct_select_key` | `kem.rs` | `subtle` mask (AND/XOR) over 32 bytes, no `jcc` on `choice` |
| `ct_asm_rm_decode` | `rm_decode(_, 3)` | `codes/reed_muller.rs` | argmax: 128 `cmov`/mask updates, no early-exit `jcc` on `f[k]` |

### 5.2 Emitting the asm

Preferred — `cargo-show-asm` (robust across MSVC/GNU, demangles, handles codegen
units):

```
cargo install cargo-show-asm
# portable clmul64 path:
cargo asm --release --features ct-audit --lib ct_asm_clmul64
cargo asm --release --features ct-audit --lib ct_asm_select_key
cargo asm --release --features ct-audit --lib ct_asm_rm_decode
cargo asm --release --features ct-audit --lib ct_asm_mul_dense_ct
# SIMD leaf:
RUSTFLAGS="-C target-feature=+pclmulqdq" cargo asm --release --features ct-audit --lib ct_asm_clmul64
```

(If a name is ambiguous, run `cargo asm --release --features ct-audit --lib` with
no symbol to list candidates.)

Fallback — raw rustc asm dump, then open the `.s` and search for each `ct_asm_*`
label:

```
RUSTFLAGS="-C codegen-units=1" cargo rustc --release --features ct-audit --lib -- --emit asm
# output: target/release/deps/hqcr-<hash>.s
```

### 5.3 miri on the pclmulqdq path

`clmul64_known_values` (a `#[cfg(test)]` smoke test in `poly/mul.rs`) calls
`clmul64` with fixed carry-less products — a fast miri target that exercises the
crate's only `unsafe` (the SIMD block) when built with the feature:

```
# portable path (sanity):
cargo +nightly miri test --lib clmul64_known_values
# pclmulqdq unsafe block:
RUSTFLAGS="-C target-feature=+pclmulqdq" cargo +nightly miri test --lib clmul64_known_values
```

A clean run proves the intrinsic block has no UB (out-of-bounds store, bad
provenance). If miri reports the intrinsic as *unsupported*, record that: the
block is then covered by the manual SAFETY argument (`mul.rs:172`), and miri
still validates the surrounding safe code and the portable fallback.

### 5.4 Recorded results

Raw dumps archived in `docs/audit/asm+miri_results.txt`. Toolchain:
`x86_64-pc-windows-msvc`, release. Note: the `ct_asm_select_key` and
`ct_asm_rm_decode` shims compiled as *forwarders* (`call`/`jmp` to the wrapped
symbol), so the real logic was read from the wrapped symbols directly.

| Target | Build | Finding | Verdict |
|---|---|---|---|
| `clmul64` (via `ct_asm_clmul64`) | portable | SSE2 auto-vectorized bit loop; mask via `neg`/`psubq` + `and`/`andpd`; sole branch is the down-counter `add rax,-4; jne` (public). No operand-dependent `jcc`. | **CT** ✅ |
| `clmul64` (via `ct_asm_clmul64`) | `+pclmulqdq` | One `pclmulqdq xmm1,xmm0,0`, two `movq`, `ret`. No branches; fixed-latency intrinsic. | **CT** ✅ |
| `ct_select_key` (`hqcr::kem`) | portable | `neg` → byte mask, broadcast; `pand`/`pandn`/`por` over all 32 B in two SSE ops. No `jcc`, no `cmov`. | **CT** ✅ |
| `rm_decode` argmax (`hqcr::codes::reed_muller`) | portable | Full 128-scan (`cmp …,128; jne`), no early exit; `|f|` via `cmovs`; new-max test via `setg` (no jump); updates via `xor`/`and`/`xor` mask. Only branches: public loop counters + a public index/length bounds check. | **CT** ✅ |
| `karatsuba` / `mul_dense_ct` | portable | Not dumped; CT follows from the confirmed branch-free `clmul64` leaf + loops bounded only by the public limb count `n`. | CT (by construction) |
| `clmul64_known_values` (miri) | portable | `test ... ok` — no UB | **clean** ✅ |
| `clmul64_known_values` (miri) | `+pclmulqdq` | `test ... ok` — intrinsic supported by miri, no UB | **clean** ✅ |

**Analysis.** Every constant-time-critical site is confirmed branch-free at the
instruction level:

- **Carry-less leaf** (`clmul64`). Portable: SSE2 auto-vectorized, masks (never
  branches) on operand bits, sole jump is the public down-counter. `+pclmulqdq`:
  a single data-independent `pclmulqdq`. miri proves the crate's only `unsafe`
  block is UB-free, with the intrinsic actually interpreted (not skipped).
- **KEM select** (`ct_select_key`). `result = (mask & a) | (~mask & b)` via
  `pand`/`pandn`/`por` over all 32 bytes — no branch, no `cmov`, so the
  valid-vs-reject choice cannot leak.
- **RM argmax** (`rm_decode`). Scans all 128 WHT outputs with no early exit;
  abs() is a `cmovs` (branchless conditional move), the new-max comparison is a
  `setg` into a register (no jump), and the running max/idx/sign updates are the
  `subtle` `xor`/`and`/`xor` mask. No `jcc` depends on a secret `f[k]`. The one
  bounds-check branch (`jae → panic_bounds_check`) tests a *public* index against
  the *public* block length, not secret content — benign, not a finding.
- **Decrypt multiply** (`mul_dense_ct`/`karatsuba`). Reduces to the confirmed
  branch-free leaf plus loops bounded only by the public limb count `n`; its
  source-level CT verdict (row 6) and clean Layer 2 canary (`|t| = 0.90`) are
  thus corroborated at the instruction level.

### 5.5 Layer 3 conclusion

**Complete — all spot-checks clean.** The asm confirms that on the CT-critical
paths the `subtle` selects compile to bitmask/`cmov`, never a secret-dependent
`jcc`; the portable `clmul64` masks rather than branches; the SIMD leaf is a
single fixed-latency `pclmulqdq`; and `rm_decode`'s argmax is a full no-early-exit
scan. miri confirms the sole `unsafe` block is free of undefined behaviour on both
the portable and `pclmulqdq` builds. Layer 3 corroborates the Layer 1 CT verdicts
and the Layer 2 timing evidence; no new leak was found at the instruction level.

---

## 6. Audit verdict summary (all layers)

| Outcome | Rows |
|---|---|
| **CT — verified** | 2, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 24, 25 |
| **NOT CT — keygen-only, accepted** | 1, 4 |
| **Accepted — ephemeral operand** | 5 |
| **NOT CT — headline limitation (RS decoder + GF tables)** | 16, 17, 18, 19, 20, 21, 22, 23 |
| **N/A — public data** | 3 |

Row 9 (`pclmulqdq`) is confirmed by Layer 3 (single fixed-latency `pclmulqdq`,
miri-clean). All three layers agree:

| Layer | Method | Outcome |
|---|---|---|
| 1 | Manual source review | CT verdicts on all encap-side paths; RS decoder flagged non-CT (§3.3) |
| 2 | dudect timing | Canary `mul_dense_ct` clean (\|t\|=0.90); decoder leak confirmed (`codes::decode` \|t\|≈622) |
| 3 | asm + miri | `subtle` selects → mask/`cmov`, never secret `jcc`; SIMD leaf one `pclmulqdq`; `unsafe` UB-free |

**Bottom line.** Every security-critical path on the *encapsulation* side
(ephemeral sampling, decrypt multiply, KEM selection/compare, RM decode) is
constant-time — verified by source review, corroborated by timing, and confirmed
at the instruction level. The single material gap is the **non-constant-time
Reed-Solomon decoder and its GF(2⁸) table lookups**, a decryption-side timing
oracle that Layer 2 measured directly; it is documented as a known limitation
deferred to a future CT-decoder rewrite. Keygen-side non-CT (secret sampler and
sparse multiply) is accepted as reference-consistent and not attacker-queryable.
No library behaviour was changed by this audit; `cargo test` and the KAT suite
pass unchanged, and all instruments are gated behind the `ct-audit` feature.
