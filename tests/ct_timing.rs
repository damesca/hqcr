//! Constant-time leakage harness — audit Layer 2 (empirical timing).
//!
//! Gated behind the `ct-audit` feature so it never compiles into a normal
//! `cargo test`. **Always run `--release`** — debug timing is meaningless for a
//! leakage test, and a debug build additionally panics inside the RS decoder on
//! random ciphertexts (`gf_inv(0)`'s `debug_assert`), which release disables.
//!
//! ```text
//! cargo test --release --features ct-audit --test ct_timing -- --nocapture
//! # tune sample count (default below):
//! CT_AUDIT_ITERS=400000 cargo test --release --features ct-audit --test ct_timing -- --nocapture
//! # measure the SIMD leaf too:
//! RUSTFLAGS="-C target-feature=+pclmulqdq" cargo test --release --features ct-audit --test ct_timing -- --nocapture
//! ```
//!
//! ## Method (dudect-style leakage detection)
//!
//! For each target we interleave two input *classes* and time one call per
//! iteration with the cycle counter (`rdtsc`, fenced with `lfence`):
//!
//!   * class **fix**  — a single fixed input (a fixed secret, or a fixed valid
//!     ciphertext);
//!   * class **rand** — a fresh random input drawn from a pre-generated pool.
//!
//! The class is chosen per iteration by a PRNG so the two are interleaved in
//! time, cancelling slow CPU-frequency / thermal drift. We then apply Welch's
//! (unequal-variance) t-test between the two classes, repeated over several
//! upper-percentile **crops** of the data — discarding the slowest X% removes
//! heavy-tailed OS-scheduling outliers that would otherwise inflate variance.
//! We report `max |t|` over the crops.
//!
//! Interpretation (per dudect): `|t| < 5` ⇒ no evidence of a leak; `|t| ≥ 10` ⇒
//! a leak signal. On a noisy laptop / Windows this is **evidence, not proof** —
//! a clean result does not certify constant-time, and a signal localizes where
//! to look. These tests therefore **log** their statistic and never fail the
//! build.
//!
//! ## Hygiene (why the numbers are meaningful)
//!
//!   * Input preparation — sampling a secret (the rejection sampler is itself
//!     variable-time!) or building a ciphertext — happens **before** the timed
//!     region, into pools. The timed closure does nothing but the operation
//!     under test, with its result `black_box`ed so it cannot be elided.
//!   * For `mul_dense_ct` the public operand `a` is held **fixed** across both
//!     classes and only the secret multiplicand `b` varies (fix vs random),
//!     isolating dependence on the secret.
//!   * **Equal-footprint pools (cache confound).** Both classes draw their
//!     operand from a pool of the *same* size, accessed at a random index — so
//!     both pay the same cache-miss cost and only the operand's *content*
//!     differs. An earlier revision reused one hot fixed operand while the rand
//!     class churned a cold pool; that address asymmetry alone made the
//!     known-CT `mul_dense_ct` canary report `|t|≈26`, which is how the confound
//!     was caught (see §4.1 of docs/audit/constant-time.md).
//!
//! ## Expected outcomes (see docs/audit/constant-time.md)
//!
//!   * `mul_dense_ct` — CT by construction (row 6); expect `|t| ≈ 0`. Serves as
//!     the **canary**: any signal here means the harness is measuring an
//!     artifact, not the operation.
//!   * `decaps` — a leak signal driven by the non-CT RS decoder (§3.3), but
//!     **diluted**: most of `decaps`'s ~7M cycles is the CT dense multiply
//!     (`u·y`) plus a near-constant re-encryption, so the decoder's
//!     data-dependence only nudges the total a fraction of a percent. Expect a
//!     weak-but-rising signal, conclusive only at large `N`.
//!   * `codes::decode` — the decoder **in isolation** (fix = clean codeword,
//!     early-return; rand = errored codeword, full BM/Chien/Forney). Here the
//!     variable decoder is most of the measured work, so the §3.3 leak shows
//!     up undiluted. Expect a **strong** signal.

#![cfg(feature = "ct-audit")]

use std::hint::black_box;
use std::time::Instant;

use hqcr::codes;
use hqcr::hash;
use hqcr::poly::mul::mul_dense_ct;
use hqcr::poly::sampling::{sample_fixed_weight, sample_uniform};
use hqcr::poly::Poly;
use hqcr::{decaps, encaps_deterministic, keygen_from_seed};
use hqcr::{Hqc128, HqcParams, SALT_BYTES, SEED_BYTES};

// ── Tunables ──────────────────────────────────────────────────────────────────

/// Default measurements per target (override with `CT_AUDIT_ITERS`).
const MUL_ITERS_DEFAULT: usize = 100_000;
const DECAPS_ITERS_DEFAULT: usize = 50_000;
/// `codes::decode` is cheap (~µs), so it can afford many more samples.
const CODES_ITERS_DEFAULT: usize = 300_000;
/// Distinct random inputs in the "rand" class pool.
const POOL: usize = 1024;
/// Throwaway iterations before measuring (warm caches / ramp frequency).
const WARMUP: usize = 2_000;
/// Upper-percentile crops for the t-test (1.0 = keep all samples).
const CROPS: [f64; 5] = [1.0, 0.999, 0.99, 0.95, 0.90];

fn iters(env_default: usize) -> usize {
    std::env::var("CT_AUDIT_ITERS")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(env_default)
}

// ── Cycle measurement ──────────────────────────────────────────────────────────

/// Time one invocation of `op` in CPU cycles. The `lfence`s keep the two
/// `rdtsc` reads from being reordered around the work. (Test-only `unsafe`; the
/// library's no-unsafe-outside-mul.rs rule is about `src/`.)
#[cfg(target_arch = "x86_64")]
#[inline]
fn measure_cycles(mut op: impl FnMut()) -> u64 {
    use core::arch::x86_64::{_mm_lfence, _rdtsc};
    unsafe {
        _mm_lfence();
        let start = _rdtsc();
        _mm_lfence();
        op();
        _mm_lfence();
        let end = _rdtsc();
        _mm_lfence();
        end.wrapping_sub(start)
    }
}

/// Portable fallback: nanoseconds instead of cycles. Noisier (lower resolution);
/// the t-test is scale-invariant so the verdict still holds, with less power.
#[cfg(not(target_arch = "x86_64"))]
#[inline]
fn measure_cycles(mut op: impl FnMut()) -> u64 {
    let start = Instant::now();
    op();
    start.elapsed().as_nanos() as u64
}

// ── Tiny PRNG (SplitMix64) — test randomness, not cryptographic ────────────────

struct SplitMix64(u64);
impl SplitMix64 {
    fn new(seed: u64) -> Self {
        Self(seed)
    }
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
    fn fill(&mut self, dst: &mut [u8]) {
        for chunk in dst.chunks_mut(8) {
            let r = self.next().to_le_bytes();
            chunk.copy_from_slice(&r[..chunk.len()]);
        }
    }
    fn seed32(&mut self) -> [u8; 32] {
        let mut s = [0u8; 32];
        self.fill(&mut s);
        s
    }
}

// ── Welch t-test + percentile cropping ─────────────────────────────────────────

fn welch_t(a: &[f64], b: &[f64]) -> f64 {
    let (na, nb) = (a.len() as f64, b.len() as f64);
    if na < 2.0 || nb < 2.0 {
        return 0.0;
    }
    let ma = a.iter().sum::<f64>() / na;
    let mb = b.iter().sum::<f64>() / nb;
    let va = a.iter().map(|x| (x - ma) * (x - ma)).sum::<f64>() / (na - 1.0);
    let vb = b.iter().map(|x| (x - mb) * (x - mb)).sum::<f64>() / (nb - 1.0);
    let denom = (va / na + vb / nb).sqrt();
    if denom == 0.0 {
        return 0.0;
    }
    (ma - mb) / denom
}

fn mean(v: &[f64]) -> f64 {
    if v.is_empty() {
        0.0
    } else {
        v.iter().sum::<f64>() / v.len() as f64
    }
}

/// Print the per-crop t-statistics and the overall verdict for one target.
fn report(name: &str, fix: &[f64], rnd: &[f64], elapsed_s: f64) {
    let mut combined: Vec<f64> = fix.iter().chain(rnd.iter()).copied().collect();
    combined.sort_by(|a, b| a.partial_cmp(b).unwrap());

    println!(
        "\n[ct-timing] {name}\n  samples: fix={} rand={}   means: fix={:.0} rand={:.0}   wall={:.1}s",
        fix.len(),
        rnd.len(),
        mean(fix),
        mean(rnd),
        elapsed_s,
    );

    let mut best = 0.0f64;
    let mut best_crop = 1.0f64;
    for &p in &CROPS {
        let idx = ((p * (combined.len() as f64 - 1.0)).round() as usize).min(combined.len() - 1);
        let thr = combined[idx];
        let fa: Vec<f64> = fix.iter().copied().filter(|&x| x <= thr).collect();
        let ra: Vec<f64> = rnd.iter().copied().filter(|&x| x <= thr).collect();
        let t = welch_t(&fa, &ra).abs();
        println!(
            "    crop {:>5.1}%  thr={:>10.0}  |t|={:>8.2}   (kept fix={} rand={})",
            p * 100.0,
            thr,
            t,
            fa.len(),
            ra.len(),
        );
        if t > best {
            best = t;
            best_crop = p;
        }
    }

    let verdict = if best >= 10.0 {
        "LEAK SIGNAL (|t| >= 10)"
    } else if best >= 5.0 {
        "INCONCLUSIVE (5 <= |t| < 10) — increase CT_AUDIT_ITERS"
    } else {
        "no leak signal (|t| < 5)"
    };
    println!(
        "  => max |t| = {:.2} at crop {:.1}%   ::   {verdict}",
        best,
        best_crop * 100.0
    );
}

// ── Target 1: mul_dense_ct (decrypt multiply u·y, secret y) ────────────────────

fn run_mul_dense_ct<P: HqcParams>(name: &str) {
    let n = iters(MUL_ITERS_DEFAULT);
    let mut prng = SplitMix64::new(0xC0FF_EE00_1234_5678);

    // Fixed public operand `a` (a dense ring element, like the ciphertext `u`).
    let a: Poly<P> = {
        let mut x = hash::xof(&[0xA1u8; 32]);
        sample_uniform::<P>(&mut x)
    };

    // Two pools of EQUAL footprint so both classes pay the same cache-miss cost
    // when the operand is fetched — only the byte *content* differs (one fixed
    // secret, replicated, vs distinct random secrets). See the Hygiene note: a
    // single hot fixed operand against a cold rand pool made this known-CT
    // canary report a spurious |t|≈26.
    let b_template: Poly<P> = {
        let mut x = hash::xof(&[0xB2u8; 32]);
        sample_fixed_weight::<P>(&mut x, P::OMEGA)
    };
    let fix_pool: Vec<Poly<P>> = (0..POOL).map(|_| b_template.clone()).collect();
    let rand_pool: Vec<Poly<P>> = (0..POOL)
        .map(|_| {
            let mut x = hash::xof(&prng.seed32());
            sample_fixed_weight::<P>(&mut x, P::OMEGA)
        })
        .collect();

    // Warm-up.
    for _ in 0..WARMUP {
        let r = mul_dense_ct::<P>(&a, &b_template);
        black_box(&r);
    }

    let mut fix = Vec::with_capacity(n / 2 + 1);
    let mut rnd = Vec::with_capacity(n / 2 + 1);

    let start = Instant::now();
    for _ in 0..n {
        let bit = prng.next();
        let class_rand = (bit & 1) == 1;
        let idx = (bit >> 1) as usize % POOL;
        let b: &Poly<P> = if class_rand { &rand_pool[idx] } else { &fix_pool[idx] };
        let c = measure_cycles(|| {
            let r = mul_dense_ct::<P>(&a, b);
            black_box(&r);
        });
        if class_rand {
            rnd.push(c as f64);
        } else {
            fix.push(c as f64);
        }
    }
    report(name, &fix, &rnd, start.elapsed().as_secs_f64());
}

// ── Target 2: decaps (exercises the non-CT RS decoder, §3.3) ───────────────────

fn run_decaps<P: HqcParams>(name: &str) {
    let n = iters(DECAPS_ITERS_DEFAULT);
    let mut prng = SplitMix64::new(0xDECA_0DEC_AB5E_1234);

    // One fixed key for the whole run (the attack model: fixed key, many cs).
    let (ek, dk) = keygen_from_seed::<P>(&[0x42u8; SEED_BYTES]);

    // class "fix": a VALID ciphertext (decode succeeds ⇒ syndromes all zero ⇒ RS
    // decoder takes its early-return path).
    let c_valid: Vec<u8> = {
        let m = vec![0x5Au8; P::K];
        let salt = [0x17u8; SALT_BYTES];
        let (_k, c) = encaps_deterministic::<P>(&ek, &m, &salt);
        c
    };
    let ct_len = c_valid.len();

    // Equal-footprint pools (see the note in run_mul_dense_ct): fix = copies of
    // the one valid ciphertext, rand = distinct random correct-length strings
    // (they pass the length check, then drive the decoder down data-dependent
    // error paths). Both fetched from same-sized cold pools at random indices,
    // so the only difference the t-test can see is the decode work itself.
    let fix_pool: Vec<Vec<u8>> = (0..POOL).map(|_| c_valid.clone()).collect();
    let rand_pool: Vec<Vec<u8>> = (0..POOL)
        .map(|_| {
            let mut c = vec![0u8; ct_len];
            prng.fill(&mut c);
            c
        })
        .collect();

    // Warm-up.
    for _ in 0..WARMUP {
        let k = decaps::<P>(&dk, &c_valid);
        black_box(&k);
    }

    let mut fix = Vec::with_capacity(n / 2 + 1);
    let mut rnd = Vec::with_capacity(n / 2 + 1);

    let start = Instant::now();
    for _ in 0..n {
        let bit = prng.next();
        let class_rand = (bit & 1) == 1;
        let idx = (bit >> 1) as usize % POOL;
        let c: &[u8] = if class_rand { &rand_pool[idx] } else { &fix_pool[idx] };
        let cyc = measure_cycles(|| {
            let k = decaps::<P>(&dk, c);
            black_box(&k);
        });
        if class_rand {
            rnd.push(cyc as f64);
        } else {
            fix.push(cyc as f64);
        }
    }
    report(name, &fix, &rnd, start.elapsed().as_secs_f64());
}

// ── Target 3: codes::decode — the RS decoder in isolation (§3.3) ───────────────
//
// Strips away the CT dense multiply and re-encryption that dilute the `decaps`
// signal, leaving the non-constant-time RMRS decoder as essentially all of the
// measured work. `fix` is a clean codeword (zero errors ⇒ RS syndromes zero ⇒
// early return); `rand` is the same codeword with a few whole RM blocks
// complemented (each is exactly one symbol error), so RS runs the full
// Berlekamp-Massey + Chien + Forney pipeline. A clean-vs-errored timing gap is
// the leak the §3.3 finding predicts.

/// Complement every bit of RM block `block` (a flat range of `N2` bits),
/// flipping that symbol's value — outside any timed region.
fn complement_block<P: HqcParams>(p: &mut Poly<P>, block: usize) {
    let start = block * P::N2;
    for i in start..start + P::N2 {
        if p.get_bit(i) == 1 {
            p.clear_bit(i);
        } else {
            p.set_bit(i);
        }
    }
}

fn run_codes_decode<P: HqcParams>(name: &str) {
    let n = iters(CODES_ITERS_DEFAULT);
    let mut prng = SplitMix64::new(0x5EED_C0DE_D0DE_1234);

    // Clean codeword for a fixed message: decodes with zero errors.
    let clean: Poly<P> = codes::encode::<P>(&vec![0x5Au8; P::K]);

    // `errs` distinct symbol errors (≤ δ ⇒ still correctable, so the full
    // decode pipeline runs to completion rather than bailing early).
    let errs = (P::DELTA / 2).max(1);

    // Equal-footprint pools: fix = copies of the clean codeword, rand = codewords
    // each carrying `errs` symbol errors at random block positions.
    let fix_pool: Vec<Poly<P>> = (0..POOL).map(|_| clean.clone()).collect();
    let rand_pool: Vec<Poly<P>> = (0..POOL)
        .map(|_| {
            let mut cw = clean.clone();
            let mut used = vec![false; P::N1];
            let mut placed = 0;
            while placed < errs {
                let b = (prng.next() as usize) % P::N1;
                if !used[b] {
                    used[b] = true;
                    complement_block::<P>(&mut cw, b);
                    placed += 1;
                }
            }
            cw
        })
        .collect();

    // Warm-up.
    for _ in 0..WARMUP {
        let r = codes::decode::<P>(&clean);
        black_box(&r);
    }

    let mut fix = Vec::with_capacity(n / 2 + 1);
    let mut rnd = Vec::with_capacity(n / 2 + 1);

    let start = Instant::now();
    for _ in 0..n {
        let bit = prng.next();
        let class_rand = (bit & 1) == 1;
        let idx = (bit >> 1) as usize % POOL;
        let cw: &Poly<P> = if class_rand { &rand_pool[idx] } else { &fix_pool[idx] };
        let cyc = measure_cycles(|| {
            let r = codes::decode::<P>(cw);
            black_box(&r);
        });
        if class_rand {
            rnd.push(cyc as f64);
        } else {
            fix.push(cyc as f64);
        }
    }
    report(name, &fix, &rnd, start.elapsed().as_secs_f64());
}

// ── Entry points ───────────────────────────────────────────────────────────────
//
// These are evidence-gathering instruments: they print a t-statistic and never
// assert, because timing on a general-purpose OS is noisy (evidence, not proof).
// Hqc128 only by default to bound runtime; add 192/256 calls if desired.

#[test]
fn ct_timing_mul_dense_ct() {
    run_mul_dense_ct::<Hqc128>("mul_dense_ct (Hqc128) — expect |t| ~ 0 (CT canary)");
}

#[test]
fn ct_timing_decaps() {
    run_decaps::<Hqc128>("decaps (Hqc128) — weak/diluted leak signal (RS decoder, §3.3)");
}

#[test]
fn ct_timing_codes_decode() {
    run_codes_decode::<Hqc128>("codes::decode (Hqc128) — RS decoder isolated, strong signal expected (§3.3)");
}
