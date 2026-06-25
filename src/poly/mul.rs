// Polynomial multiplication in R = F2[X]/(X^n - 1) — the hot path.
//
// Two modes:
//
//   mul_sparse_dense (Mode A): one operand has low Hamming weight (ω or ωr).
//     Iterates over the set-bit positions of the sparse operand and
//     XOR-adds a cyclic rotation of the dense operand into the accumulator.
//     Cost: O(weight · N_WORDS).
//     Used in: keygen (h·y, h·r2), encrypt (h·r2, s·r2).
//     No CT requirement — positions come from the public XOF.
//
//   mul_dense_ct (Mode B): both operands are arbitrary, but the second
//     operand's bit positions must not be leaked via timing (it may be a
//     secret key component y).
//     Uses the same rotate-and-accumulate strategy but extracts bits with
//     a branchless mask so every word is visited regardless of its value.
//     Cost: O(N · N_WORDS) — much slower, but only called in Decrypt.
//     CT requirement: no data-dependent branches on `b`'s bits.
//
// After either multiplication call reduce() on the result to fold
// the overflow bits back into [0, N).
//
// ── Optimization layers ───────────────────────────────────────────────────────
//
// Mode A (mul_sparse_dense) stays at the O(weight·N_WORDS) rotate-and-accumulate
// baseline — Karatsuba cannot beat a linear-in-N algorithm, so the layers below
// apply only to the O(N²) Mode B dense×dense path (the decrypt hot path).
//
// L0 [DONE] — portable word-level baseline.
//   Mode A (this file) and the reduce_wide / shift_xor_wide infrastructure.
//
// L1 [DONE] — Karatsuba over GF(2)[X], always compiled and always used by
//   mul_dense_ct. Implemented at the LIMB level (split the operand word arrays
//   on a limb boundary h = ⌈n/2⌉) so every sub-operand stays word-aligned and
//   recombination needs no sub-word shifts. The classic 3-multiply identity
//   A·B = z0 + (z1 − z0 − z2)·βʰ + z2·β²ʰ over GF(2) (subtraction = XOR),
//   bottoming out in a schoolbook base case of carry-less word multiplies.
//   See `karatsuba` below.
//
// L2 [DONE] — SIMD carry-less multiply via pclmulqdq (x86-64), the leaf word
//   multiply `clmul64`. Selected at COMPILE TIME by
//   `#[cfg(target_feature = "pclmulqdq")]`: build with
//   `RUSTFLAGS="-C target-feature=+pclmulqdq"` (or `-C target-cpu=native`) to
//   enable it. A branch-free portable `clmul64` is always present as the
//   fallback on every other target. The pclmulqdq intrinsic block is the ONLY
//   `unsafe` in the crate. Reference: the official HQC AVX2 C implementation
//   (pqc-hqc.org) uses _mm_clmulepi64_si128 in its vect_mul for the same job.

use super::{Poly, MAX_N_WORDS};
use crate::params::HqcParams;

// ── Wide accumulator ──────────────────────────────────────────────────────────
//
// Both multiply modes accumulate shifted copies of one operand into a
// double-width, NON-wrapping buffer, then reduce that buffer modulo X^N - 1 in
// a single fold. This is the correct way to multiply in R:
//
//   * A product of two degree-(<N) polynomials has degree ≤ 2N-2, so it needs
//     up to 2N bits — a buffer of 2·N_WORDS words always suffices (plus a small
//     slack so the fold's `acc[word+1]` read never runs past the end).
//   * Shifting WITHOUT a `% N_WORDS` wrap keeps every bit at its true integer
//     position p+rot. (The previous implementation wrapped word indices modulo
//     N_WORDS, which reduces modulo X^(N_WORDS·64)-1 ≠ X^N-1 and silently
//     misplaced any bit with p+rot ≥ N_WORDS·64 by N_WORDS·64 - N positions.)
//   * Reduction mod X^N-1 means X^N ≡ 1, so bit (t+N) folds onto bit t. Since
//     all bits live in [0, 2N), each low bit t∈[0,N) receives exactly one fold
//     contribution: result_bit[t] = acc_bit[t] XOR acc_bit[t+N].
//
// reduce_wide / shift_xor_wide are shared: Mode A uses them directly; Mode B's
// Karatsuba writes its 2N-bit product into the same wide layout, then reduces.

/// Width of the non-wrapping product accumulator, in u64 words.
/// 2·MAX_N_WORDS covers the ≤2N-bit product; +2 gives slack for the fold read.
const WIDE_WORDS: usize = 2 * MAX_N_WORDS + 2;

/// XOR `src` (an N-bit ring element in its low `N_WORDS` words) shifted left by
/// `rot` bit positions into the wide accumulator `acc`. No modular wraparound:
/// bits land at their true positions [rot, rot+N).
#[inline(never)] // single call target for the inner loop's branch predictor
fn shift_xor_wide<P: HqcParams>(acc: &mut [u64; WIDE_WORDS], src: &[u64; MAX_N_WORDS], rot: usize) {
    let nw = P::N_WORDS;
    let word_shift = rot >> 6; // rot / 64
    let bit_shift = rot & 63; // rot % 64

    if bit_shift == 0 {
        for i in 0..nw {
            acc[i + word_shift] ^= src[i];
        }
    } else {
        let right_shift = 64 - bit_shift;
        for i in 0..nw {
            acc[i + word_shift] ^= src[i] << bit_shift;
            acc[i + word_shift + 1] ^= src[i] >> right_shift;
        }
    }
    // Max index touched: (nw-1) + word_shift + 1. With rot < N, word_shift ≤
    // (N-1)/64 ≤ nw-1, so the index is ≤ 2·nw-1 < WIDE_WORDS. No bounds issue.
}

/// Fold the wide accumulator (bits in [0, 2N)) into a reduced `Poly<P>` modulo
/// X^N - 1: `result_bit[t] = acc_bit[t] XOR acc_bit[t+N]` for t in [0, N).
fn reduce_wide<P: HqcParams>(acc: &[u64; WIDE_WORDS]) -> Poly<P> {
    let nw = P::N_WORDS;
    let n = P::N;
    let word = n / 64; // word holding bit N
    let off = n % 64; // bit offset of N within that word

    let mut out = Poly::<P>::zero();
    for w in 0..nw {
        // Extract the 64 bits of `acc` starting at bit position N + w·64,
        // i.e. at word (word + w), bit offset `off`.
        let folded = if off == 0 {
            acc[word + w]
        } else {
            (acc[word + w] >> off) | (acc[word + w + 1] << (64 - off))
        };
        out.words[w] = acc[w] ^ folded;
    }

    // Clear any bits at/above N in the top word (the fold can leave junk there).
    if off != 0 {
        out.words[nw - 1] &= (1u64 << off) - 1;
    }
    out
}

// ── Mode A: sparse × dense (non-CT) ──────────────────────────────────────────

/// Multiply `sparse` (low weight) by `dense` in R.
///
/// Collects the set-bit positions of `sparse` via a regular branch (safe
/// since sparse positions are public — drawn from a public XOF).
/// For each set position `p`, XOR-shifts `dense` by `p` into the wide
/// accumulator, then reduces once.
///
/// Callers: keygen and encrypt (h·y, h·r2, s·r2).
pub fn mul_sparse_dense<P: HqcParams>(sparse: &Poly<P>, dense: &Poly<P>) -> Poly<P> {
    let mut acc = [0u64; WIDE_WORDS];

    for word_idx in 0..P::N_WORDS {
        let mut word = sparse.words[word_idx];
        while word != 0 {
            // Extract lowest set bit position.
            let lsb = word.trailing_zeros() as usize;
            let pos = word_idx * 64 + lsb;
            if pos < P::N {
                shift_xor_wide::<P>(&mut acc, &dense.words, pos);
            }
            word &= word - 1; // clear lowest set bit
        }
    }

    reduce_wide::<P>(&acc)
}

// ── GF(2)[X] limb multiplication (Karatsuba + carry-less leaf) ────────────────
//
// Constant-time: a full product visits every limb regardless of operand
// values — there is no data-dependent branch or memory access anywhere in
// `karatsuba` or in either `clmul64`. So Mode B stays constant-time with
// respect to both operands (in particular the secret `y` in Decrypt).
// pclmulqdq is itself a data-independent-latency instruction.

/// Carry-less 64×64 → 128-bit multiply over GF(2)[X]. Returns `(lo, hi)`, the
/// low and high 64 bits of the product. pclmulqdq path — one instruction.
#[cfg(all(target_arch = "x86_64", target_feature = "pclmulqdq"))]
#[inline]
fn clmul64(a: u64, b: u64) -> (u64, u64) {
    use core::arch::x86_64::{__m128i, _mm_clmulepi64_si128, _mm_set_epi64x, _mm_storeu_si128};
    // SAFETY: gated on `target_feature = "pclmulqdq"`, so the carry-less
    // multiply intrinsic is guaranteed available; `_mm_set_epi64x` /
    // `_mm_storeu_si128` are x86-64 SSE2 baseline. The store targets the local
    // 16-byte `out` array, so no out-of-bounds access is possible.
    unsafe {
        let xa = _mm_set_epi64x(0, a as i64);
        let xb = _mm_set_epi64x(0, b as i64);
        let prod = _mm_clmulepi64_si128(xa, xb, 0x00);
        let mut out = [0u64; 2];
        _mm_storeu_si128(out.as_mut_ptr() as *mut __m128i, prod);
        (out[0], out[1])
    }
}

/// Carry-less 64×64 → 128-bit multiply over GF(2)[X], portable branch-free
/// fallback. Constant-time: iterates all 64 bits with a mask, never branching
/// on operand values (only on the public loop counter).
#[cfg(not(all(target_arch = "x86_64", target_feature = "pclmulqdq")))]
#[inline]
fn clmul64(a: u64, b: u64) -> (u64, u64) {
    let mut lo = 0u64;
    let mut hi = 0u64;
    let mut i = 0;
    while i < 64 {
        // mask = all-ones iff bit i of b is set (no branch on the bit value).
        let mask = ((b >> i) & 1).wrapping_neg();
        lo ^= (a << i) & mask;
        if i != 0 {
            // High word receives the bits of (a << i) that spill past bit 63;
            // i ∈ 1..64 keeps the shift in range.
            hi ^= (a >> (64 - i)) & mask;
        }
        i += 1;
    }
    (lo, hi)
}

/// Limb count at or below which `karatsuba` switches to its schoolbook base
/// case. Tunable knob — measure with `cargo bench --bench bench -- poly_mul`.
const KARATSUBA_THRESHOLD: usize = 16;

/// Scratch-arena size (u64 words) for `karatsuba`. The recursion needs
/// ≈ 4·(n + log₂ n) words in total (4·⌈n/2⌉ locally per level plus the
/// geometric tail of the sub-calls); this bound sits comfortably above the
/// worst case at `MAX_N_WORDS`.
const KARATSUBA_SCRATCH_WORDS: usize = 4 * MAX_N_WORDS + 256;

/// Multiply the `n`-limb GF(2)[X] polynomials `a` and `b`, **overwriting**
/// `out[0..2n]` with the `2n`-limb product. `scratch` is a caller-provided
/// workspace; `KARATSUBA_SCRATCH_WORDS` words always suffice.
///
/// Limb-level Karatsuba split on `h = ⌈n/2⌉` (high half `n − h ≤ h`):
///   z0 = a_lo·b_lo,  z2 = a_hi·b_hi,
///   mid = (a_lo+a_hi)·(b_lo+b_hi) − z0 − z2     (− is XOR over GF(2)),
///   product = z0 + mid·βʰ + z2·β²ʰ              (β = 2⁶⁴),
/// applied recursively down to a schoolbook base case of `clmul64`s.
fn karatsuba(out: &mut [u64], a: &[u64], b: &[u64], n: usize, scratch: &mut [u64]) {
    if n <= KARATSUBA_THRESHOLD {
        // Schoolbook: out[0..2n] = Σ_{i,j} clmul64(a[i], b[j]) placed at limb i+j.
        for slot in out[..2 * n].iter_mut() {
            *slot = 0;
        }
        for i in 0..n {
            let ai = a[i];
            for j in 0..n {
                let (lo, hi) = clmul64(ai, b[j]);
                out[i + j] ^= lo;
                out[i + j + 1] ^= hi;
            }
        }
        return;
    }

    let h = n.div_ceil(2); // low-half limb count (ceil); high half = n - h ≤ h
    let nhi = n - h;

    let (a_lo, a_hi) = (&a[0..h], &a[h..n]);
    let (b_lo, b_hi) = (&b[0..h], &b[h..n]);

    // Partition the arena: sum_a (h) | sum_b (h) | mid (2h) | sub (recursion).
    let (sum_a, rest) = scratch.split_at_mut(h);
    let (sum_b, rest) = rest.split_at_mut(h);
    let (mid, sub) = rest.split_at_mut(2 * h);

    // sum_a = a_lo + a_hi, sum_b = b_lo + b_hi (high halves shorter ⇒ XOR low end).
    sum_a.copy_from_slice(a_lo);
    for i in 0..nhi {
        sum_a[i] ^= a_hi[i];
    }
    sum_b.copy_from_slice(b_lo);
    for i in 0..nhi {
        sum_b[i] ^= b_hi[i];
    }

    // z0 → out[0..2h], z2 → out[2h..2n] (recursive calls overwrite their range).
    karatsuba(&mut out[0..2 * h], a_lo, b_lo, h, sub);
    karatsuba(&mut out[2 * h..2 * n], a_hi, b_hi, nhi, sub);

    // mid = sum_a·sum_b, then subtract (XOR) z0 and z2 already sitting in `out`.
    karatsuba(mid, sum_a, sum_b, h, sub);
    for i in 0..2 * h {
        mid[i] ^= out[i]; // − z0
    }
    for i in 0..2 * nhi {
        mid[i] ^= out[2 * h + i]; // − z2
    }

    // product += mid·βʰ (overlaps the high half of z0 and the low half of z2).
    for i in 0..2 * h {
        out[h + i] ^= mid[i];
    }
}

// ── Mode B: dense × dense, constant-time on `b` ──────────────────────────────

/// Multiply `a` by `b` in R, constant-time with respect to both operands'
/// bit patterns. Used in Decrypt where `b = y` is the secret key.
///
/// Computes the full `2N`-bit GF(2)[X] product with limb-level Karatsuba
/// (`karatsuba`, leaf `clmul64`) into the wide accumulator, then folds it
/// modulo `X^N − 1` with `reduce_wide`. Because a full product touches every
/// limb unconditionally, there is no data-dependent control flow — timing is
/// independent of the secret. Cost: O(N^log₂3) ≈ O(N^1.585), the Karatsuba
/// improvement over the O(N²) bit-level baseline.
///
/// Relies on the crate-wide `Poly` invariant that bits ≥ N are zero in both
/// operands (held by every sampler and ring op), so the symmetric full-limb
/// product equals `a · b` in R.
pub fn mul_dense_ct<P: HqcParams>(a: &Poly<P>, b: &Poly<P>) -> Poly<P> {
    let nw = P::N_WORDS;

    // Debug-only guard on the invariant the equivalence above depends on.
    debug_assert!(
        {
            let top = P::N & 63;
            if top == 0 {
                true
            } else {
                let junk = !((1u64 << top) - 1);
                (a.words[nw - 1] & junk) == 0 && (b.words[nw - 1] & junk) == 0
            }
        },
        "mul_dense_ct operands must have zero bits at/above N"
    );

    let mut acc = [0u64; WIDE_WORDS];
    let mut scratch = [0u64; KARATSUBA_SCRATCH_WORDS];
    karatsuba(
        &mut acc[..2 * nw],
        &a.words[..nw],
        &b.words[..nw],
        nw,
        &mut scratch,
    );
    reduce_wide::<P>(&acc)
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::{Hqc128, Hqc192, Hqc256};

    // Helper: build a Poly with exactly the given bit positions set.
    fn from_positions<P: HqcParams>(positions: &[usize]) -> Poly<P> {
        let mut p = Poly::<P>::zero();
        for &i in positions {
            p.set_bit(i);
        }
        p
    }

    // ── Basic algebraic identities ────────────────────────────────────────────

    #[test]
    fn multiply_by_zero_is_zero() {
        let mut a = Poly::<Hqc128>::zero();
        a.set_bit(0);
        a.set_bit(100);
        let zero = Poly::<Hqc128>::zero();
        let r = mul_sparse_dense::<Hqc128>(&zero, &a);
        assert_eq!(r.hamming_weight(), 0);
        let r2 = mul_sparse_dense::<Hqc128>(&a, &zero);
        assert_eq!(r2.hamming_weight(), 0);
    }

    #[test]
    fn multiply_by_one_is_identity() {
        // "1" in R is the polynomial with only the constant term set (bit 0).
        let one = from_positions::<Hqc128>(&[0]);
        let mut a = Poly::<Hqc128>::zero();
        a.set_bit(5);
        a.set_bit(1000);
        a.set_bit(17000);
        let r = mul_sparse_dense::<Hqc128>(&one, &a);
        assert_eq!(r, a);
        let r2 = mul_sparse_dense::<Hqc128>(&a, &one);
        assert_eq!(r2, a);
    }

    #[test]
    fn multiply_by_x_is_cyclic_shift() {
        // Multiplying by X shifts all exponents up by 1 mod N.
        // bit i → bit (i+1) mod N.
        let x = from_positions::<Hqc128>(&[1]); // X^1
        let a = from_positions::<Hqc128>(&[0, 100, Hqc128::N - 1]);
        let r = mul_sparse_dense::<Hqc128>(&x, &a);
        // bit 0 → bit 1
        // bit 100 → bit 101
        // bit N-1 → bit 0 (cyclic wrap)
        assert_eq!(r.get_bit(1), 1);
        assert_eq!(r.get_bit(101), 1);
        assert_eq!(r.get_bit(0), 1, "X^(N-1) * X should wrap to bit 0");
        assert_eq!(r.hamming_weight(), 3);
    }

    // ── Absolute correctness vs naive cyclic convolution ─────────────────────
    //
    // This is the test the old code lacked: it compares against an independent
    // O(weight²) reference, exercising LARGE rotations (positions near N) that
    // trigger the wraparound path. Commutativity/distributivity could not catch
    // the previous bug because they compared the buggy multiply against itself.

    /// Naive product in R: bit (s+d) mod N toggled for each pair of set bits.
    fn naive_mul<P: HqcParams>(a_pos: &[usize], b_pos: &[usize]) -> Poly<P> {
        let mut acc = vec![0u8; P::N];
        for &s in a_pos {
            for &d in b_pos {
                let k = (s + d) % P::N;
                acc[k] ^= 1;
            }
        }
        let mut p = Poly::<P>::zero();
        for (i, &bit) in acc.iter().enumerate() {
            if bit == 1 {
                p.set_bit(i);
            }
        }
        p
    }

    fn check_against_naive<P: HqcParams>(a_pos: &[usize], b_pos: &[usize]) {
        let a = from_positions::<P>(a_pos);
        let b = from_positions::<P>(b_pos);
        let expected = naive_mul::<P>(a_pos, b_pos);
        assert_eq!(mul_sparse_dense::<P>(&a, &b), expected, "Mode A vs naive");
        assert_eq!(mul_dense_ct::<P>(&b, &a), expected, "Mode B vs naive");
    }

    #[test]
    fn matches_naive_large_rotations_128() {
        // Positions near N force pos + (N-1) ≫ N_WORDS·64 — the wraparound path.
        check_against_naive::<Hqc128>(
            &[0, 1, 12_345, Hqc128::N - 1, Hqc128::N - 2],
            &[0, 7, 17_000, Hqc128::N - 1],
        );
    }

    #[test]
    fn matches_naive_large_rotations_192() {
        check_against_naive::<Hqc192>(
            &[0, 3, 30_000, Hqc192::N - 1],
            &[5, 200, 35_000, Hqc192::N - 5],
        );
    }

    #[test]
    fn matches_naive_large_rotations_256() {
        check_against_naive::<Hqc256>(
            &[0, 9, 50_000, Hqc256::N - 1],
            &[1, 64, 57_000, Hqc256::N - 3],
        );
    }

    #[test]
    fn matches_naive_single_high_bit_128() {
        // X^(N-1) · X^(N-1) = X^(2N-2) = X^(N-2): a minimal wraparound check.
        check_against_naive::<Hqc128>(&[Hqc128::N - 1], &[Hqc128::N - 1]);
    }

    // ── Commutativity ─────────────────────────────────────────────────────────

    #[test]
    fn sparse_dense_commutativity_hqc128() {
        use crate::poly::sampling::sample_fixed_weight;
        use sha3::{
            digest::{ExtendableOutput, Update},
            Shake256,
        };

        for seed in 0u8..10 {
            let mut xof = {
                let mut h = Shake256::default();
                h.update(&[seed, 0]);
                h.finalize_xof()
            };
            let a = sample_fixed_weight::<Hqc128>(&mut xof, Hqc128::OMEGA);
            let mut xof2 = {
                let mut h = Shake256::default();
                h.update(&[seed, 1]);
                h.finalize_xof()
            };
            let b = sample_fixed_weight::<Hqc128>(&mut xof2, Hqc128::OMEGA_R);

            let ab = mul_sparse_dense::<Hqc128>(&a, &b);
            let ba = mul_sparse_dense::<Hqc128>(&b, &a);
            assert_eq!(ab, ba, "commutativity failed for seed={seed}");
        }
    }

    #[test]
    fn sparse_dense_commutativity_hqc256() {
        use crate::poly::sampling::sample_fixed_weight;
        use sha3::{
            digest::{ExtendableOutput, Update},
            Shake256,
        };

        let mut xof = {
            let mut h = Shake256::default();
            h.update(b"hqc256-seed");
            h.finalize_xof()
        };
        let a = sample_fixed_weight::<Hqc256>(&mut xof, Hqc256::OMEGA);
        let mut xof2 = {
            let mut h = Shake256::default();
            h.update(b"hqc256-seed2");
            h.finalize_xof()
        };
        let b = sample_fixed_weight::<Hqc256>(&mut xof2, Hqc256::OMEGA_R);

        let ab = mul_sparse_dense::<Hqc256>(&a, &b);
        let ba = mul_sparse_dense::<Hqc256>(&b, &a);
        assert_eq!(ab, ba);
    }

    // ── Mode A vs Mode B agreement ────────────────────────────────────────────

    #[test]
    fn sparse_dense_matches_dense_ct() {
        // Both modes must produce identical results given the same inputs.
        use crate::poly::sampling::{sample_fixed_weight, sample_uniform};
        use sha3::{
            digest::{ExtendableOutput, Update},
            Shake256,
        };

        for seed in 0u8..5 {
            let mut xof = {
                let mut h = Shake256::default();
                h.update(&[seed]);
                h.finalize_xof()
            };
            // sparse operand (secret-like)
            let sparse = sample_fixed_weight::<Hqc128>(&mut xof, Hqc128::OMEGA);
            // dense operand (ciphertext component)
            let dense = sample_uniform::<Hqc128>(&mut xof);

            let r_a = mul_sparse_dense::<Hqc128>(&sparse, &dense);
            let r_b = mul_dense_ct::<Hqc128>(&dense, &sparse);
            assert_eq!(r_a, r_b, "Mode A vs Mode B mismatch for seed={seed}");
        }
    }

    // ── Karatsuba vs the original bit-level multiply ──────────────────────────
    //
    // The L0 dense×dense reference (the implementation mul_dense_ct replaced):
    // bit-by-bit, branchless on `b`. Kept here as an independent oracle so the
    // Karatsuba/clmul path is checked against it across random inputs. A
    // divergence here is exactly the kind of bug that would break KAT.

    fn mul_dense_ct_bitwise<P: HqcParams>(a: &Poly<P>, b: &Poly<P>) -> Poly<P> {
        let mut acc = [0u64; WIDE_WORDS];
        let nw = P::N_WORDS;
        for pos in 0..P::N {
            let word_idx = pos >> 6;
            let bit_idx = pos & 63;
            let bit = (b.words[word_idx] >> bit_idx) & 1;
            let mask = bit.wrapping_neg();
            let word_shift = pos >> 6;
            let bit_shift = pos & 63;
            if bit_shift == 0 {
                for i in 0..nw {
                    acc[i + word_shift] ^= a.words[i] & mask;
                }
            } else {
                let right_shift = 64 - bit_shift;
                for i in 0..nw {
                    acc[i + word_shift] ^= (a.words[i] << bit_shift) & mask;
                    acc[i + word_shift + 1] ^= (a.words[i] >> right_shift) & mask;
                }
            }
        }
        reduce_wide::<P>(&acc)
    }

    fn karatsuba_matches_bitwise<P: HqcParams>(tag: &[u8]) {
        use crate::poly::sampling::sample_uniform;
        use sha3::{
            digest::{ExtendableOutput, Update},
            Shake256,
        };

        for seed in 0u8..8 {
            let mut xof = {
                let mut h = Shake256::default();
                h.update(tag);
                h.update(&[seed]);
                h.finalize_xof()
            };
            // Two dense, uniformly random operands (tail bits ≥ N masked off by
            // sample_uniform, so the invariant mul_dense_ct relies on holds).
            let a = sample_uniform::<P>(&mut xof);
            let b = sample_uniform::<P>(&mut xof);

            let fast = mul_dense_ct::<P>(&a, &b);
            let reference = mul_dense_ct_bitwise::<P>(&a, &b);
            assert_eq!(
                fast, reference,
                "Karatsuba vs bit-level mismatch, seed={seed}"
            );

            // Commutativity of the full product, for good measure.
            let fast_rev = mul_dense_ct::<P>(&b, &a);
            assert_eq!(fast, fast_rev, "dense product not commutative, seed={seed}");
        }
    }

    #[test]
    fn karatsuba_matches_bitwise_128() {
        karatsuba_matches_bitwise::<Hqc128>(b"k128");
    }
    #[test]
    fn karatsuba_matches_bitwise_192() {
        karatsuba_matches_bitwise::<Hqc192>(b"k192");
    }
    #[test]
    fn karatsuba_matches_bitwise_256() {
        karatsuba_matches_bitwise::<Hqc256>(b"k256");
    }

    // ── Distributivity ────────────────────────────────────────────────────────

    #[test]
    fn distributivity_add_mul() {
        // a * (b + c) == a*b + a*c
        let a = from_positions::<Hqc128>(&[0, 1, 5]);
        let b = from_positions::<Hqc128>(&[2, 7]);
        let c = from_positions::<Hqc128>(&[3, 8, 100]);

        let bc = b.add(&c);
        let lhs = mul_sparse_dense::<Hqc128>(&a, &bc);
        let ab = mul_sparse_dense::<Hqc128>(&a, &b);
        let ac = mul_sparse_dense::<Hqc128>(&a, &c);
        let rhs = ab.add(&ac);
        assert_eq!(lhs, rhs);
    }

    // ── Reduce correctness ────────────────────────────────────────────────────

    #[test]
    fn result_has_no_overflow_bits_after_mul() {
        use crate::poly::sampling::sample_fixed_weight;
        use sha3::{
            digest::{ExtendableOutput, Update},
            Shake256,
        };

        let mut xof = {
            let mut h = Shake256::default();
            h.update(b"overflow-check");
            h.finalize_xof()
        };
        let a = sample_fixed_weight::<Hqc128>(&mut xof, Hqc128::OMEGA);
        let b = sample_fixed_weight::<Hqc128>(&mut xof, Hqc128::OMEGA_R);
        let r = mul_sparse_dense::<Hqc128>(&a, &b);

        let last_bit = Hqc128::N & 63;
        let mask = (1u64 << last_bit) - 1;
        assert_eq!(
            r.words[Hqc128::N_WORDS - 1] & !mask,
            0,
            "overflow bits not cleared by reduce()"
        );
    }

    // ── Parameter set compilation checks ─────────────────────────────────────

    #[test]
    fn all_three_param_sets_mul() {
        let a128 = from_positions::<Hqc128>(&[0, 1]);
        let b128 = from_positions::<Hqc128>(&[0, 2]);
        let _ = mul_sparse_dense::<Hqc128>(&a128, &b128);

        let a192 = from_positions::<Hqc192>(&[0, 1]);
        let b192 = from_positions::<Hqc192>(&[0, 2]);
        let _ = mul_sparse_dense::<Hqc192>(&a192, &b192);

        let a256 = from_positions::<Hqc256>(&[0, 1]);
        let b256 = from_positions::<Hqc256>(&[0, 2]);
        let _ = mul_sparse_dense::<Hqc256>(&a256, &b256);
    }

    // ── clmul64 known values (miri target for the pclmulqdq unsafe block) ──────

    #[test]
    fn clmul64_known_values() {
        // Carry-less (GF(2)[X]) products — independent of which `clmul64` impl is
        // selected, so the same test exercises the portable mask loop and, under
        // `+pclmulqdq`, the `_mm_clmulepi64_si128` intrinsic. Small + fast, the
        // ideal `cargo miri test` target for the crate's only `unsafe`.
        assert_eq!(clmul64(0, 0), (0, 0));
        assert_eq!(clmul64(1, 0xDEAD_BEEF), (0xDEAD_BEEF, 0)); // multiply by 1
        assert_eq!(clmul64(2, 2), (4, 0)); // x · x = x^2
        assert_eq!(clmul64(0b11, 0b11), (0b101, 0)); // (x+1)^2 = x^2 + 1 over GF(2)
        assert_eq!(clmul64(1 << 63, 2), (0, 1)); // x^63 · x = x^64 → hi bit 0
        assert_eq!(clmul64(u64::MAX, 1), (u64::MAX, 0));
    }
}
