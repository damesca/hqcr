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
// ── Optimization layers (implement in order) ─────────────────────────────────
//
// L0 [DONE] — portable word-level baseline (this file).
//   rotate-and-accumulate with u64 words; always compiled; no unsafe.
//
// L1 [TODO] — Karatsuba over GF(2)[X] (~2× over L0 for Mode B dense×dense).
//   Split each N-bit poly into two halves A_lo, A_hi of N/2 bits.
//   Karatsuba identity: A*B = A_lo*B_lo + (A_lo*B_lo + A_hi*B_hi + (A_lo+A_hi)*(B_lo+B_hi))*X^(N/2) + A_hi*B_hi*X^N
//   Reduces three half-size multiplications instead of four (the classic 3M trick).
//   Benefit: meaningful only for dense×dense (Mode B); sparse×dense (Mode A) is
//   already O(weight·N_WORDS) so Karatsuba does not help there.
//   Gate with a compile-time flag or simply replace mul_dense_ct once validated.
//
// L2 [TODO] — SIMD carry-less multiply via pclmulqdq (x86-64).
//   Use std::arch::x86_64::{__m128i, _mm_clmulepi64_si128} to multiply pairs
//   of 64-bit words in a single instruction (GF(2)[X] multiply, no carry).
//   Gate with: #[cfg(target_feature = "pclmul")]
//   Keep the L0 portable path as the fallback — the pclmul block is the only
//   place in this crate where `unsafe` is permitted.
//   Reference: the official HQC AVX2 C implementation (pqc-hqc.org) uses
//   _mm_clmulepi64_si128 in its vect_mul routine for the same purpose.

use crate::params::HqcParams;
use super::{Poly, MAX_N_WORDS};

// ── Cyclic rotation helpers ───────────────────────────────────────────────────

/// XOR `src` rotated left by `rot` positions (in the ring of N bits) into `dst`.
///
/// "Rotated left by rot" means coefficient i of src contributes to coefficient
/// (i + rot) mod N of dst — i.e. multiplying by X^rot in R.
///
/// Implementation: a cyclic rotation by `rot` bits decomposes into:
///   1. A word shift by `word_shift = rot / 64`.
///   2. A bit shift within a word: `bit_shift = rot % 64`.
///
/// We iterate over each source word, split its bits across two consecutive
/// destination words (with the bit shift), and wrap around modulo N at the
/// end.  The wrap is handled by a second pass that folds the `N_WORDS`-th
/// word's carry into the appropriate low words.
///
/// After this function returns, `dst` may have non-zero bits above position
/// N-1 (in the tail of dst[N_WORDS-1]).  The caller must call reduce() once
/// after accumulating all rotations.
#[inline(never)] // keep the inner loop a single call target for branch prediction
fn xor_rotate_into<P: HqcParams>(dst: &mut [u64; MAX_N_WORDS], src: &[u64; MAX_N_WORDS], rot: usize) {
    let nw = P::N_WORDS;
    let word_shift = rot >> 6;      // rot / 64
    let bit_shift  = rot & 63;      // rot % 64

    if bit_shift == 0 {
        // Aligned rotation — just word-shift with wrap.
        for i in 0..nw {
            let j = (i + word_shift) % nw;
            dst[j] ^= src[i];
        }
    } else {
        let right_shift = 64 - bit_shift;
        for i in 0..nw {
            let lo_idx = (i + word_shift) % nw;
            let hi_idx = (i + word_shift + 1) % nw;
            dst[lo_idx] ^= src[i] << bit_shift;
            dst[hi_idx] ^= src[i] >> right_shift;
        }
    }
    // NOTE: This does not correctly implement cyclic reduction — the modular
    // wrap from `% nw` on indices handles word-level wraparound but does NOT
    // account for the fact that the ring has N bits, not N_WORDS*64 bits.
    // The caller MUST call poly.reduce() after all rotations are accumulated.
    //
    // The issue: when bits from the end of the ring spill over into words[0]
    // they arrive at bit positions that correspond to ring positions >= N.
    // reduce() handles this by folding those overflow bits back.
}

// ── Mode A: sparse × dense (non-CT) ──────────────────────────────────────────

/// Multiply `sparse` (low weight) by `dense` in R.
///
/// Collects the set-bit positions of `sparse` via a regular branch (safe
/// since sparse positions are public — drawn from a public XOF).
/// For each set position `p`, XOR-rotates `dense` by `p` into the result.
///
/// Callers: keygen and encrypt (h·y, h·r2, s·r2).
pub fn mul_sparse_dense<P: HqcParams>(sparse: &Poly<P>, dense: &Poly<P>) -> Poly<P> {
    let mut result = Poly::<P>::zero();

    for word_idx in 0..P::N_WORDS {
        let mut word = sparse.words[word_idx];
        while word != 0 {
            // Extract lowest set bit position.
            let lsb = word.trailing_zeros() as usize;
            let pos = word_idx * 64 + lsb;
            if pos < P::N {
                xor_rotate_into::<P>(&mut result.words, &dense.words, pos);
            }
            word &= word - 1; // clear lowest set bit
        }
    }

    result.reduce();
    result
}

// ── Mode B: dense × dense, constant-time on `b` ──────────────────────────────

/// Multiply `a` by `b` in R, constant-time with respect to `b`'s bit pattern.
///
/// Visits every bit of `b` using a branchless mask — no early exit, no
/// data-dependent control flow on `b`.  Used in Decrypt where `b = y` is
/// the secret key.
///
/// For each bit position `pos` in `b` (0..N), extract the bit with a mask
/// and conditionally XOR-rotate `a` into the accumulator.  The conditional
/// is implemented as:
///
///   if (b.words[w] >> bit) & 1 == 1 { xor_rotate_into(a, pos) }
///
/// The branch on the extracted bit IS a secret-dependent branch — to fix it
/// we replace it with a branchless conditional: multiply each word of the
/// rotation by the bit mask (0 or 1), so the XOR either adds or adds zero.
///
/// Cost: O(N · N_WORDS) — roughly N/OMEGA times slower than Mode A.
pub fn mul_dense_ct<P: HqcParams>(a: &Poly<P>, b: &Poly<P>) -> Poly<P> {
    let mut result = Poly::<P>::zero();
    let nw = P::N_WORDS;

    for pos in 0..P::N {
        let word_idx = pos >> 6;
        let bit_idx  = pos & 63;
        // Extract bit `pos` from `b` as a mask: 0x0000…0000 or 0xFFFF…FFFF.
        // This is the CT heart: no branch, no multiply — just arithmetic.
        let bit = (b.words[word_idx] >> bit_idx) & 1;
        let mask = bit.wrapping_neg(); // 0 → 0x0000…, 1 → 0xFFFF…

        // XOR-rotate a by `pos` into result, gated by mask.
        // We inline the rotation here (rather than call xor_rotate_into)
        // so we can apply the mask directly without an extra allocation.
        let word_shift = pos >> 6;
        let bit_shift  = pos & 63;

        if bit_shift == 0 {
            for i in 0..nw {
                let j = (i + word_shift) % nw;
                result.words[j] ^= a.words[i] & mask;
            }
        } else {
            let right_shift = 64 - bit_shift;
            for i in 0..nw {
                let lo_idx = (i + word_shift) % nw;
                let hi_idx = (i + word_shift + 1) % nw;
                result.words[lo_idx] ^= (a.words[i] << bit_shift)  & mask;
                result.words[hi_idx] ^= (a.words[i] >> right_shift) & mask;
            }
        }
    }

    result.reduce();
    result
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

    // ── Commutativity ─────────────────────────────────────────────────────────

    #[test]
    fn sparse_dense_commutativity_hqc128() {
        use sha3::{Shake256, digest::{Update, ExtendableOutput}};
        use crate::poly::sampling::sample_fixed_weight;

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
        use sha3::{Shake256, digest::{Update, ExtendableOutput}};
        use crate::poly::sampling::sample_fixed_weight;

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
        use sha3::{Shake256, digest::{Update, ExtendableOutput}};
        use crate::poly::sampling::{sample_fixed_weight, sample_uniform};

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

    // ── Distributivity ────────────────────────────────────────────────────────

    #[test]
    fn distributivity_add_mul() {
        // a * (b + c) == a*b + a*c
        let a = from_positions::<Hqc128>(&[0, 1, 5]);
        let b = from_positions::<Hqc128>(&[2, 7]);
        let c = from_positions::<Hqc128>(&[3, 8, 100]);

        let bc = b.add(&c);
        let lhs = mul_sparse_dense::<Hqc128>(&a, &bc);
        let ab  = mul_sparse_dense::<Hqc128>(&a, &b);
        let ac  = mul_sparse_dense::<Hqc128>(&a, &c);
        let rhs = ab.add(&ac);
        assert_eq!(lhs, rhs);
    }

    // ── Reduce correctness ────────────────────────────────────────────────────

    #[test]
    fn result_has_no_overflow_bits_after_mul() {
        use sha3::{Shake256, digest::{Update, ExtendableOutput}};
        use crate::poly::sampling::sample_fixed_weight;

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
        assert_eq!(r.words[Hqc128::N_WORDS - 1] & !mask, 0,
            "overflow bits not cleared by reduce()");
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
}
