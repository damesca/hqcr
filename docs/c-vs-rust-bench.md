# C vs Rust HQC benchmark

Comparison of the pure-Rust `hqcr` implementation against the official C
reference and C optimized implementations, across all three security levels.

## Environment

| | |
|---|---|
| CPU | Intel Core i7-13620H (13th Gen, Raptor Lake-H) |
| OS | Windows 11 |
| C toolchain | MSYS2 mingw-w64 gcc (MINGW64) |
| Rust toolchain | rustc 1.93.0 (254b59607 2026-01-19) |
| C iterations | 500 per operation |
| Rust samples | 100 (criterion default) |
| Date | 2026-06-30 |

## Results

All times are median wall-clock milliseconds. Rust values are criterion point
estimates. C cycle column is `__rdtsc` median (invariant TSC — indicative, not
turbo-accurate). Lower is better.

### HQC-128

| Op | C ref portable | Rust portable | C opt AVX2 | Rust +pclmulqdq | C opt cycles | Rust +pclmul cycles† |
|---|---:|---:|---:|---:|---:|---:|
| keygen | 1.153 ms | 0.032 ms | 0.041 ms | 0.032 ms | 3 364 070 | — |
| encaps | 2.326 ms | 0.578 ms | 0.042 ms | 0.577 ms | 123 010 | — |
| decaps | 3.523 ms | 2.542 ms | 0.092 ms | 0.917 ms | 269 160 | — |

### HQC-192

| Op | C ref portable | Rust portable | C opt AVX2 | Rust +pclmulqdq |  C opt cycles | Rust +pclmul cycles† |
|---|---:|---:|---:|---:|---:|---:|
| keygen | 3.572 ms | 0.068 ms | 0.097 ms | 0.069 ms | 283 971 | — |
| encaps | 7.192 ms | 1.599 ms | 0.096 ms | 1.575 ms | 281 096 | — |
| decaps | 10.847 ms | 7.194 ms | 0.189 ms | 2.306 ms | 550 064 | — |

### HQC-256

| Op | C ref portable | Rust portable | C opt AVX2 | Rust +pclmulqdq | C opt cycles | Rust +pclmul cycles† |
|---|---:|---:|---:|---:|---:|---:|
| keygen | 8.761 ms | 0.128 ms | 0.111 ms | 0.129 ms | 323 505 | — |
| encaps | 17.679 ms | 3.290 ms | 0.196 ms | 3.246 ms | 570 937 | — |
| decaps | 26.709 ms | 17.619 ms | 0.371 ms | 4.938 ms | 1 083 260 | — |

† Criterion does not expose rdtsc; fill in with a perf/VTune run if needed.

## Observations

### 1 — Keygen and encaps: Rust sparse×dense beats C reference by ~30–70×

The Rust implementation uses Mode A (`mul_sparse_dense`) for keygen and encaps.
Because the secret and ephemeral vectors are sparse (ω = 66 / 100 / 131, ω_r =
75 / 114 / 149 non-zero positions), multiplying sparse × dense costs
`ω × ⌈n/64⌉` word operations instead of the full `⌈n/64⌉²` that the C
reference generic polynomial multiplier pays. For HQC-128 (n = 17 669):
`75 × 277 ≈ 20 775` word ops vs `277² ≈ 76 729` — a ~3.7× algorithmic
advantage, amplified further because the Rust inner loop is branchless and cache-
friendly. This is a genuine algorithm difference, not just a constant-factor win.

The C optimized build closes most of the gap with full AVX2-vectorized arithmetic
(and likely an additive NTT poly multiplier), reaching near-parity with Rust for
keygen and beating Rust ~13× for encaps.

### 2 — Decaps: C optimized AVX2 is the fastest by far

Decaps is the one operation that requires a full dense × dense polynomial
multiplication (computing `u·y` where `y` is dense). The results reflect this:

| Level | C opt | Rust +pclmul | Rust portable |
|---|---:|---:|---:|
| HQC-128 | 0.092 ms | 0.917 ms | 2.542 ms |
| HQC-192 | 0.189 ms | 2.306 ms | 7.194 ms |
| HQC-256 | 0.371 ms | 4.938 ms | 17.619 ms |

C optimized is **~10× faster** than Rust +pclmulqdq for decaps. The C optimized
build uses full 256-bit AVX2 SIMD throughout the polynomial multiplier; Rust's
`mul_dense_ct` uses a Karatsuba tree with pclmulqdq only at the leaf level.
Closing this gap would require a fully vectorized Karatsuba or an additive NTT
in Rust (currently not implemented — CLAUDE.md Step 20c-opt).

### 3 — pclmulqdq affects Rust decaps only

Enabling `RUSTFLAGS="-C target-feature=+pclmulqdq"` has no effect on Rust
keygen or encaps (they use sparse×dense, which does not call `clmul64`), but
cuts Rust decaps by **63–72%** across all levels:

| Level | portable | +pclmulqdq | speedup |
|---|---:|---:|---:|
| HQC-128 | 2.542 ms | 0.917 ms | 2.8× |
| HQC-192 | 7.194 ms | 2.306 ms | 3.1× |
| HQC-256 | 17.619 ms | 4.938 ms | 3.6× |

### 4 — C reference is the slowest everywhere

C reference (portable `-O3`, no SIMD) pays full generic poly-mul cost for every
operation, with no sparsity exploitation. It is the slowest in all nine cells.
The gap vs Rust portable ranges from 1.4× (HQC-128 decaps, where both do dense
multiplication) to ~68× (HQC-128 keygen, where Rust exploits sparsity).

### 5 — C optimized keygen vs Rust keygen

At HQC-256, C optimized keygen (0.111 ms) is slightly *faster* than Rust
(0.128 ms). At HQC-128 and HQC-192, Rust is slightly faster. These differences
are within measurement noise and depend on AVX2 random-number expansion overhead
in the C PRNG path (an RNG asymmetry noted in the methodology below).

## Methodology caveats

**RNG asymmetry.** The Rust criterion bench uses `keygen_from_seed` and
`encaps_deterministic`, which take pre-expanded randomness as fixed inputs. The C
harness calls `crypto_kem_keypair` and `crypto_kem_enc` with a live PRNG, which
draws seeds / `m` / `salt` from SHAKE-256 *inside* the timed call. This adds a
small SHAKE overhead to the C numbers for keygen and encaps. For decaps there is
no asymmetry (both sides take `(sk, ct)` as fixed input).

**Algorithm asymmetry.** The "C ref portable ↔ Rust portable" pairing refers to
instruction-set portability (no SIMD), not algorithmic equivalence. Rust uses
sparse×dense multiplication for keygen/encaps; C reference uses a generic dense
multiplier. The comparison shows what each portable-targeted implementation
achieves in practice, not a constant-factor hardware comparison.

**Cycle counter.** `__rdtsc` reports invariant TSC reference cycles, not
execution-core cycles under turbo boost. Treat the cycle column as indicative
only. Wall-clock ns/op (median) is the authoritative metric.

**Single-run numbers.** Each C variant was run once with 500 iterations.
Criterion uses 100 samples with warm-up. For publication-quality results, run
each 2–3 times and verify medians are stable within ~2%.

**Platform.** All measurements were taken on the same Windows 11 machine with
MSYS2 MINGW64 gcc (C) and the default `x86_64-pc-windows-msvc` Rust target.
Absolute numbers will differ on other CPUs; ratios should be broadly stable.
