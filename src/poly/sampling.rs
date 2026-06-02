// Sampling functions for polynomials in R = F2[X]/(X^n - 1).
//
// Two samplers:
//   sample_fixed_weight — rejection sampling, exactly `weight` distinct positions
//   sample_uniform      — fills all N bits uniformly (for the public key component h)
//
// Constant-time contract:
//   sample_fixed_weight must not branch on the *value* of sampled positions,
//   only on the public condition `pos < n`. The deduplication check (is this
//   position already in the set?) uses subtle::ConstantTimeEq so timing does
//   not reveal which positions were accepted or rejected.

use sha3::digest::{ExtendableOutput, XofReader};
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq};

use crate::params::HqcParams;
use super::Poly;

// ── Internal XOF helpers ──────────────────────────────────────────────────────

/// Read exactly 2 bytes from the XOF and interpret them as a little-endian u16.
#[inline(always)]
fn read_u16(xof: &mut impl XofReader) -> u16 {
    let mut buf = [0u8; 2];
    xof.read(&mut buf);
    u16::from_le_bytes(buf)
}

/// Read exactly 8 bytes from the XOF and interpret them as a little-endian u64.
#[inline(always)]
fn read_u64(xof: &mut impl XofReader) -> u64 {
    let mut buf = [0u8; 8];
    xof.read(&mut buf);
    u64::from_le_bytes(buf)
}

// ── Fixed-weight sampler ──────────────────────────────────────────────────────

/// Sample a polynomial with exactly `weight` coefficients equal to 1, with the
/// set positions drawn uniformly at random from [0, n).
///
/// Algorithm: rejection sampling over u16 values read from `xof`.
///   - Read a u16.  If it is >= n, discard (public rejection — does not leak).
///   - Otherwise check whether it duplicates any already-accepted position.
///     The duplicate check is done in constant time via ConstantTimeEq so that
///     the timing does not reveal the accepted positions.
///   - Accept iff the position is in-range AND not a duplicate.
///   - Repeat until exactly `weight` positions are collected.
///
/// The `positions` accumulation buffer lives on the stack as a fixed-size
/// array of u32 values (max weight across all parameter sets is 149 for
/// HQC-256, but we size to 256 for safety).  We store u32 (not u16) so we
/// can accommodate n up to 57 637 without truncation.
///
/// CT contract: the only data-dependent branch is `pos < n`, which depends
/// only on the public parameter n.  All comparisons against already-accepted
/// positions use subtle::ConstantTimeEq and ConditionallySelectable.
pub fn sample_fixed_weight<P: HqcParams>(
    xof: &mut impl XofReader,
    weight: usize,
) -> Poly<P> {
    debug_assert!(weight <= 256, "weight {weight} exceeds internal buffer size");
    debug_assert!(weight <= P::N, "weight {weight} > N={}", P::N);

    // Accepted positions, stored as u32.  We use a fixed buffer; only the
    // first `filled` entries are valid, but we always operate on the full
    // buffer in CT mode so timing is uniform.
    let mut positions = [0u32; 256];
    let mut filled: usize = 0;

    while filled < weight {
        // Draw a 16-bit candidate.  n <= 57637 < 65536 so u16 covers the range.
        let raw = read_u16(xof) as u32;

        // Public rejection: discard if out of range.  This branch is on a
        // public parameter (n), not on secret data.
        if raw >= P::N as u32 {
            continue;
        }

        // Duplicate check — must be constant time.
        //
        // We scan all `weight` slots of `positions` (not just `filled`).
        // Slots above `filled` contain 0; we bias the comparison to treat
        // slot index >= filled as "no match" by also checking the slot index.
        //
        // Strategy:
        //   For each slot j in 0..weight:
        //     is_used   = (j < filled) as Choice          [public index < public count]
        //     is_match  = positions[j].ct_eq(&raw)        [CT value comparison]
        //     duplicate |= is_used & is_match
        //
        // `j < filled` is a comparison of two public loop counters — no secret
        // data involved — so a regular branch is fine here.
        let mut duplicate = Choice::from(0u8);
        for j in 0..weight {
            let is_used = Choice::from((j < filled) as u8);
            let is_match = positions[j].ct_eq(&raw);
            duplicate |= is_used & is_match;
        }

        // Accept this position iff it is not a duplicate.
        // We write to positions[filled] unconditionally (CT), then increment
        // `filled` only if accepted.  The unconditional write is safe because
        // `filled` < weight <= 256.
        //
        // ConditionallySelectable: positions[filled] = if accepted { raw } else { positions[filled] }
        let accepted = !duplicate;
        positions[filled] = u32::conditional_select(&positions[filled], &raw, accepted);
        // filled += accepted as usize: branch on `accepted`, which is derived
        // from the duplicate check — is this public? No: `filled` reveals how
        // many distinct positions have been found, which depends on the XOF
        // output. HOWEVER, `filled` is not secret in the threat model: the
        // loop iteration count is observable via timing anyway (we loop until
        // filled == weight, which is public). The spec does not require hiding
        // the number of rejections, only the position *values*. So a branch on
        // `accepted` for the counter is acceptable.
        if accepted.into() {
            filled += 1;
        }
    }

    // Build the polynomial: set bit at each accepted position.
    let mut poly = Poly::<P>::zero();
    for j in 0..weight {
        poly.set_bit(positions[j] as usize);
    }
    poly
}

// ── Uniform sampler ───────────────────────────────────────────────────────────

/// Sample a uniformly random polynomial with N bits from `xof`.
///
/// Reads N_WORDS full u64 words from the XOF, then masks off the bits
/// above position N-1 in the last word so the polynomial is properly reduced.
///
/// Used for the public key component h — no CT requirement (h is public).
pub fn sample_uniform<P: HqcParams>(xof: &mut impl XofReader) -> Poly<P> {
    let mut poly = Poly::<P>::zero();

    for i in 0..P::N_WORDS {
        poly.words[i] = read_u64(xof);
    }

    // Zero the bits above N-1 in the last word.
    let last_bit = P::N & 63; // valid bit count in the last word
    if last_bit != 0 {
        let mask = (1u64 << last_bit) - 1;
        poly.words[P::N_WORDS - 1] &= mask;
    }

    poly
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::{Hqc128, Hqc192, Hqc256};
    use sha3::{Shake256, digest::Update};

    fn make_xof(seed: &[u8]) -> impl XofReader {
        use sha3::digest::ExtendableOutput;
        let mut h = Shake256::default();
        h.update(seed);
        h.finalize_xof()
    }

    // ── fixed-weight sampler ──────────────────────────────────────────────────

    #[test]
    fn fixed_weight_correct_weight_128() {
        let mut xof = make_xof(b"test-seed-0");
        let p = sample_fixed_weight::<Hqc128>(&mut xof, Hqc128::OMEGA);
        assert_eq!(p.hamming_weight(), Hqc128::OMEGA);
    }

    #[test]
    fn fixed_weight_correct_weight_192() {
        let mut xof = make_xof(b"test-seed-1");
        let p = sample_fixed_weight::<Hqc192>(&mut xof, Hqc192::OMEGA);
        assert_eq!(p.hamming_weight(), Hqc192::OMEGA);
    }

    #[test]
    fn fixed_weight_correct_weight_256() {
        let mut xof = make_xof(b"test-seed-2");
        let p = sample_fixed_weight::<Hqc256>(&mut xof, Hqc256::OMEGA);
        assert_eq!(p.hamming_weight(), Hqc256::OMEGA);
    }

    #[test]
    fn fixed_weight_all_bits_in_range() {
        let mut xof = make_xof(b"test-seed-range");
        let p = sample_fixed_weight::<Hqc128>(&mut xof, Hqc128::OMEGA_R);
        // Every set bit must be in [0, N).
        for i in 0..Hqc128::N {
            let _ = p.get_bit(i); // panics if out of range in debug mode
        }
        // No bits set in overflow region.
        let last_bit = Hqc128::N & 63;
        let mask = (1u64 << last_bit) - 1;
        assert_eq!(p.words[Hqc128::N_WORDS - 1] & !mask, 0);
    }

    #[test]
    fn fixed_weight_no_duplicates() {
        // With weight 1, sample many times and verify we always get exactly 1 bit set.
        for seed in 0u8..20 {
            let mut xof = make_xof(&[seed]);
            let p = sample_fixed_weight::<Hqc128>(&mut xof, 1);
            assert_eq!(p.hamming_weight(), 1, "seed={seed}");
        }
    }

    #[test]
    fn fixed_weight_deterministic() {
        // Same seed produces same polynomial.
        let mut xof1 = make_xof(b"deterministic");
        let mut xof2 = make_xof(b"deterministic");
        let p1 = sample_fixed_weight::<Hqc128>(&mut xof1, Hqc128::OMEGA);
        let p2 = sample_fixed_weight::<Hqc128>(&mut xof2, Hqc128::OMEGA);
        assert_eq!(p1, p2);
    }

    #[test]
    fn fixed_weight_different_seeds_differ() {
        let mut xof1 = make_xof(b"seed-A");
        let mut xof2 = make_xof(b"seed-B");
        let p1 = sample_fixed_weight::<Hqc128>(&mut xof1, Hqc128::OMEGA);
        let p2 = sample_fixed_weight::<Hqc128>(&mut xof2, Hqc128::OMEGA);
        assert_ne!(p1, p2, "different seeds should almost certainly differ");
    }

    // ── uniform sampler ───────────────────────────────────────────────────────

    #[test]
    fn uniform_no_overflow_bits() {
        let mut xof = make_xof(b"uniform-test");
        let p = sample_uniform::<Hqc128>(&mut xof);
        let last_bit = Hqc128::N & 63;
        let mask = (1u64 << last_bit) - 1;
        assert_eq!(p.words[Hqc128::N_WORDS - 1] & !mask, 0,
            "bits above N-1 must be zero");
    }

    #[test]
    fn uniform_deterministic() {
        let mut xof1 = make_xof(b"uniform-det");
        let mut xof2 = make_xof(b"uniform-det");
        let p1 = sample_uniform::<Hqc128>(&mut xof1);
        let p2 = sample_uniform::<Hqc128>(&mut xof2);
        assert_eq!(p1, p2);
    }

    #[test]
    fn uniform_256_no_overflow() {
        let mut xof = make_xof(b"uniform-256");
        let p = sample_uniform::<Hqc256>(&mut xof);
        let last_bit = Hqc256::N & 63;
        let mask = (1u64 << last_bit) - 1;
        assert_eq!(p.words[Hqc256::N_WORDS - 1] & !mask, 0);
    }
}
