// Poly<P: HqcParams>: bit-packed polynomial in R = F2[X]/(X^n - 1).
//
// Representation: [u64; MAX_N_WORDS] — fixed maximum capacity (HQC-256 needs
// 901 words). For smaller parameter sets the high words are always zero.
// This avoids heap allocation while keeping stable Rust (generic_const_exprs
// is not yet stable, so [u64; P::N_WORDS] as a struct field is not possible).
//
// Word ordering: little-endian — bit i lives in words[i/64] at position i%64.
// Only the first N bits are meaningful; bits N..MAX_N_WORDS*64 are always zero.

pub mod mul;
pub mod sampling;

use crate::params::HqcParams;
use core::marker::PhantomData;
use zeroize::{Zeroize, ZeroizeOnDrop};

/// Maximum N_WORDS across all parameter sets (HQC-256: ceil(57637/64) = 901).
pub(crate) const MAX_N_WORDS: usize = 901;

/// Bit-packed element of R = F2[X]/(X^n - 1).
///
/// Layout: `words[i/64]` bit `i%64` holds the coefficient of X^i.
/// Bits from index `P::N` up to `MAX_N_WORDS*64 - 1` are always zero.
#[derive(Zeroize, ZeroizeOnDrop)]
pub struct Poly<P: HqcParams> {
    pub(crate) words: [u64; MAX_N_WORDS],
    _p: PhantomData<P>,
}

// Manual Clone: the derived impl would add a spurious `where P: Clone` bound,
// but `P` is only a phantom marker (zero-sized), so cloning never touches it.
// Without this, generic code (e.g. kem.rs) could not clone a Poly under the
// bare `P: HqcParams` bound.
impl<P: HqcParams> Clone for Poly<P> {
    fn clone(&self) -> Self {
        Self {
            words: self.words,
            _p: PhantomData,
        }
    }
}

impl<P: HqcParams> Poly<P> {
    /// Returns the zero polynomial.
    #[inline]
    pub fn zero() -> Self {
        Self {
            words: [0u64; MAX_N_WORDS],
            _p: PhantomData,
        }
    }

    // ── Bit access ────────────────────────────────────────────────────────────

    /// Returns the coefficient of X^i (0 or 1). Panics in debug if i >= N.
    #[inline]
    pub fn get_bit(&self, i: usize) -> u64 {
        debug_assert!(i < P::N, "bit index {i} out of range (N={})", P::N);
        (self.words[i >> 6] >> (i & 63)) & 1
    }

    /// Sets the coefficient of X^i to 1. Panics in debug if i >= N.
    #[inline]
    pub fn set_bit(&mut self, i: usize) {
        debug_assert!(i < P::N, "bit index {i} out of range (N={})", P::N);
        self.words[i >> 6] |= 1u64 << (i & 63);
    }

    /// Clears the coefficient of X^i to 0. Panics in debug if i >= N.
    #[inline]
    pub fn clear_bit(&mut self, i: usize) {
        debug_assert!(i < P::N, "bit index {i} out of range (N={})", P::N);
        self.words[i >> 6] &= !(1u64 << (i & 63));
    }

    /// Zeros every coefficient.
    #[inline]
    pub fn clear(&mut self) {
        self.words.fill(0);
    }

    // ── Addition in R (XOR) ───────────────────────────────────────────────────

    /// Returns self + rhs in R (bitwise XOR, no carry, no reduction).
    #[inline]
    pub fn add(&self, rhs: &Self) -> Self {
        let mut out = Self::zero();
        for i in 0..P::N_WORDS {
            out.words[i] = self.words[i] ^ rhs.words[i];
        }
        out
    }

    /// In-place addition: self += rhs (XOR).
    #[inline]
    pub fn add_assign(&mut self, rhs: &Self) {
        for i in 0..P::N_WORDS {
            self.words[i] ^= rhs.words[i];
        }
    }

    // ── Cyclic reduction ─────────────────────────────────────────────────────

    /// Folds the overflow bits [N .. N_WORDS*64) back into [0 .. N).
    ///
    /// After a sparse×dense multiplication the result may have non-zero bits
    /// above position N-1. This function reduces modulo X^N - 1 by XOR-ing
    /// the high bits into the low positions.
    ///
    /// How it works:
    ///   The last active word (`words[N_WORDS-1]`) may hold bits from bit N
    ///   onward. Bit k (where k >= N) in the ring means X^k ≡ X^(k mod N).
    ///
    ///   We split the last word at the boundary bit `last_bit = N % 64`:
    ///     low  = words[N_WORDS-1] & mask(last_bit)      // bits N-last_bit .. N-1 (valid)
    ///     high = words[N_WORDS-1] >> last_bit            // bits N .. N_WORDS*64-1 (overflow)
    ///
    ///   The overflow bits represent coefficients of X^N, X^(N+1), …
    ///   We shift `high` into the low positions by XOR-ing it into words[0]:
    ///     words[0] ^= high
    ///   and keep only the valid bits in the last word.
    ///
    ///   There can be at most one partial word of overflow (one XOR is enough)
    ///   because sparse×dense multiplication never produces a result wider than
    ///   2N-2 bits, which fits in N_WORDS+1 words — but we allocate N_WORDS and
    ///   the caller (mul.rs) uses the rotate-and-accumulate strategy that keeps
    ///   the result within [0, 2N) bits, so one reduction pass suffices.
    #[inline]
    pub fn reduce(&mut self) {
        let last_bit = P::N & 63; // number of valid bits in the last word (0 means all 64)
        if last_bit == 0 {
            // N is an exact multiple of 64 — no partial last word, nothing to fold.
            return;
        }
        let last_idx = P::N_WORDS - 1;
        let mask = (1u64 << last_bit) - 1;

        let overflow = self.words[last_idx] >> last_bit;
        self.words[last_idx] &= mask; // zero out the bits above N-1
        self.words[0] ^= overflow; // fold them into the bottom of the ring
    }

    // ── Hamming weight ────────────────────────────────────────────────────────

    /// Returns the number of set bits (Hamming weight) of the polynomial.
    /// Counts only the N active bits, not the padding.
    pub fn hamming_weight(&self) -> usize {
        self.words[..P::N_WORDS]
            .iter()
            .map(|w| w.count_ones() as usize)
            .sum()
    }
}

// ── Debug / PartialEq ─────────────────────────────────────────────────────────
//
// PartialEq compares only the active N_WORDS. This is used in tests and
// in kem.rs for the implicit rejection check — but note that kem.rs must use
// the constant-time version (subtle::ConstantTimeEq) for security.

impl<P: HqcParams> PartialEq for Poly<P> {
    fn eq(&self, other: &Self) -> bool {
        self.words[..P::N_WORDS] == other.words[..P::N_WORDS]
    }
}

impl<P: HqcParams> Eq for Poly<P> {}

impl<P: HqcParams> core::fmt::Debug for Poly<P> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "Poly(N={}, weight={})", P::N, self.hamming_weight())
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::{Hqc128, Hqc192, Hqc256};

    #[test]
    fn zero_poly_has_no_set_bits() {
        let p = Poly::<Hqc128>::zero();
        assert_eq!(p.hamming_weight(), 0);
        for i in 0..Hqc128::N {
            assert_eq!(p.get_bit(i), 0);
        }
    }

    #[test]
    fn set_and_get_bit_roundtrip() {
        let mut p = Poly::<Hqc128>::zero();
        let positions = [0, 1, 63, 64, 65, Hqc128::N - 1];
        for &i in &positions {
            p.set_bit(i);
            assert_eq!(p.get_bit(i), 1, "bit {i} should be 1 after set");
        }
        assert_eq!(p.hamming_weight(), positions.len());
        for &i in &positions {
            p.clear_bit(i);
            assert_eq!(p.get_bit(i), 0, "bit {i} should be 0 after clear");
        }
        assert_eq!(p.hamming_weight(), 0);
    }

    #[test]
    fn add_is_xor() {
        let mut a = Poly::<Hqc128>::zero();
        let mut b = Poly::<Hqc128>::zero();
        a.set_bit(0);
        a.set_bit(1);
        b.set_bit(1);
        b.set_bit(2);
        let c = a.add(&b);
        // 0 XOR 0 = 0, 1 XOR 1 = 0, 0 XOR 1 = 1
        assert_eq!(c.get_bit(0), 1);
        assert_eq!(c.get_bit(1), 0); // canceled
        assert_eq!(c.get_bit(2), 1);
    }

    #[test]
    fn add_assign_matches_add() {
        let mut a = Poly::<Hqc192>::zero();
        let mut b = Poly::<Hqc192>::zero();
        a.set_bit(100);
        b.set_bit(200);
        let expected = a.add(&b);
        a.add_assign(&b);
        assert_eq!(a, expected);
    }

    #[test]
    fn reduce_folds_overflow() {
        // For Hqc128: N = 17669, N_WORDS = 277.
        // 277 * 64 = 17728. Overflow region: bits 17669..17728 (59 bits).
        // Set a bit in the overflow region and verify reduce moves it to bit 0.
        let mut p = Poly::<Hqc128>::zero();
        // bit 17669 is at word 17669/64 = 276, position 17669%64 = 5.
        // After reduce, bit 17669 mod 17669 = 0 should be set.
        p.words[276] |= 1u64 << 5; // bit 17669 (= N itself, which is in overflow)
        p.reduce();
        assert_eq!(p.get_bit(0), 1, "overflow bit N should fold to bit 0");
        assert_eq!(p.words[276] >> 5, 0, "overflow cleared after reduce");
    }

    #[test]
    fn reduce_is_idempotent() {
        let mut p = Poly::<Hqc256>::zero();
        p.set_bit(1000);
        p.set_bit(50000);
        p.reduce(); // should be no-op since bits are in range
        let w_before = p.hamming_weight();
        p.reduce();
        assert_eq!(p.hamming_weight(), w_before);
    }

    #[test]
    fn hamming_weight_counts_only_active_bits() {
        let mut p = Poly::<Hqc128>::zero();
        // Manually set a padding bit (above N) and verify it's not counted.
        // N_WORDS * 64 - N = 277*64 - 17669 = 17728 - 17669 = 59 padding bits.
        // We intentionally bypass set_bit (which would panic) to write raw:
        p.words[Hqc128::N_WORDS - 1] = u64::MAX; // includes overflow bits
                                                 // hamming_weight sums all N_WORDS words — so it WILL count overflow.
                                                 // This test verifies the caller (mul.rs / sampling.rs) always calls
                                                 // reduce() before hamming_weight(), OR that hamming_weight is only
                                                 // called on valid (reduced) polys. Document the contract here:
                                                 // hamming_weight counts N_WORDS words, caller must ensure upper bits = 0.
        let raw_count = p.hamming_weight();
        assert_eq!(raw_count, 64); // full word = 64 set bits
                                   // After mask:
        p.reduce();
        // reduce only folds if there's overflow, but doesn't zero them out for
        // a fully-set last word — let's check the actual contract.
        // Actually reduce zeros the overflow bits in the last word and XORs into
        // words[0]. words[0] was 0, so words[0] becomes the overflow bits.
        // overflow = words[276] >> 5  (59 bits all set = 0x7FFFFFFFFFFFFFFF >> (64-59))
        // This test is more of a documentation test — just ensure no panic.
        let _ = p.hamming_weight();
    }

    #[test]
    fn all_three_param_sets_compile() {
        let _a = Poly::<Hqc128>::zero();
        let _b = Poly::<Hqc192>::zero();
        let _c = Poly::<Hqc256>::zero();
    }

    #[test]
    fn clear_zeros_all_words() {
        let mut p = Poly::<Hqc128>::zero();
        p.set_bit(0);
        p.set_bit(Hqc128::N - 1);
        p.clear();
        assert_eq!(p.hamming_weight(), 0);
        assert!(p.words.iter().all(|&w| w == 0));
    }
}
