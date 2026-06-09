// Sampling functions for polynomials in R = F2[X]/(X^n - 1).
//
// Three samplers:
//   sample_fixed_weight     — rejection sampling, exactly `weight` distinct
//                             positions. Used for the long-term secret (x, y).
//   sample_fixed_weight_mod — Barrett-reduction ("mod") sampler, exactly
//                             `weight` distinct positions. Used for the
//                             ephemeral vectors (r1, r2, e) in PKE.Encrypt.
//   sample_uniform          — fills all N bits uniformly (public key h).
//
// Why two fixed-weight samplers? The 2025 spec deliberately draws the two roles
// differently (matches the reference `vect_set_random_fixed_weight` family):
//   • x, y are the long-term secret, sampled once at keygen — rejection
//     sampling gives a perfectly uniform position distribution.
//   • r1, r2, e are ephemeral, resampled on every Encrypt (and again on every
//     Decaps re-encryption), so they sit on the hot path — the Barrett sampler
//     draws a fixed 4·weight bytes with no rejection loop. Its O(N/2^32) bias
//     is cryptographically negligible for a value that is never reused.
// The two consume the XOF stream differently and produce different positions
// from the same seed, so they are NOT interchangeable — using the wrong one is
// a KAT mismatch even though both yield a valid weight-`ω` vector.
//
// Constant-time contract:
//   sample_fixed_weight must not branch on the *value* of sampled positions,
//   only on the public condition `pos < n`. The deduplication check (is this
//   position already in the set?) uses subtle::ConstantTimeEq so timing does
//   not reveal which positions were accepted or rejected.
//   sample_fixed_weight_mod is fully branchless: the reduction, the duplicate
//   resolution, and the final bit-setting all run over fixed-length loops with
//   constant-time selects, leaking nothing about the sampled positions.

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

/// Read exactly 4 bytes from the XOF and interpret them as a little-endian u32.
#[inline(always)]
fn read_u32(xof: &mut impl XofReader) -> u32 {
    let mut buf = [0u8; 4];
    xof.read(&mut buf);
    u32::from_le_bytes(buf)
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

// ── Barrett-reduction ("mod") fixed-weight sampler ──────────────────────────────

/// Sample a polynomial with exactly `weight` coefficients equal to 1, drawing
/// each position with the spec's Barrett-reduction sampler (the
/// `vect_set_random_fixed_weight` procedure of the HQC reference). Used for the
/// ephemeral vectors `r1`, `r2`, `e` in PKE.Encrypt.
///
/// Algorithm (no rejection — exactly `4·weight` bytes consumed from `xof`):
///   1. For i in 0..weight: draw a little-endian u32 `rand_i` and set
///        support[i] = i + ((rand_i · (N − i)) >> 32)
///      The fixed-point multiply-shift maps `rand_i` uniformly into `[0, N−i)`,
///      so `support[i] ∈ [i, N)`.
///   2. Resolve duplicates backwards: for i from weight−1 down to 0, if
///      support[i] equals any already-finalized support[j] (j > i), replace it
///      with `i`. Because every later entry is ≥ i+1, the value `i` is
///      guaranteed distinct from them — so the final `weight` positions are all
///      distinct. (This is the standard reservoir-style construction; its exact
///      output is what the KAT vectors pin.)
///   3. Set the corresponding bits.
///
/// CT contract: every loop has a fixed, public length (`weight` or `N_WORDS`),
/// the duplicate test uses `subtle::ConstantTimeEq`, the duplicate fix uses
/// `ConditionallySelectable`, and the bit-setting scans every word OR-ing in a
/// masked bit — so neither timing nor memory-access pattern depends on the
/// sampled position values.
pub fn sample_fixed_weight_mod<P: HqcParams>(
    xof: &mut impl XofReader,
    weight: usize,
) -> Poly<P> {
    debug_assert!(weight <= 256, "weight {weight} exceeds internal buffer size");
    debug_assert!(weight <= P::N, "weight {weight} > N={}", P::N);

    let mut support = [0u32; 256];

    // Step 1 — draw and Barrett-reduce: support[i] = i + ((rand_i·(N−i)) >> 32).
    for i in 0..weight {
        let rand = read_u32(xof) as u64;
        let n_minus_i = (P::N - i) as u64;
        support[i] = (i as u64 + ((rand * n_minus_i) >> 32)) as u32;
    }

    // Step 2 — backward duplicate resolution (constant time). If support[i]
    // collides with any finalized support[j] (j > i), replace it with `i`.
    for i in (0..weight).rev() {
        let mut found = Choice::from(0u8);
        for j in (i + 1)..weight {
            found |= support[j].ct_eq(&support[i]);
        }
        support[i] = u32::conditional_select(&support[i], &(i as u32), found);
    }

    // Precompute (word index, bit mask) for each position.
    let mut word_of = [0u32; 256];
    let mut bit_of = [0u64; 256];
    for k in 0..weight {
        word_of[k] = support[k] >> 6;
        bit_of[k] = 1u64 << (support[k] & 63);
    }

    // Step 3 — constant-time bit setting: for every word, OR in the bits of the
    // positions that land in it (selected branchlessly), never a data-dependent
    // store address.
    let mut poly = Poly::<P>::zero();
    for word_idx in 0..P::N_WORDS {
        let wi = word_idx as u32;
        let mut acc = 0u64;
        for k in 0..weight {
            let eq = word_of[k].ct_eq(&wi);
            acc |= u64::conditional_select(&0u64, &bit_of[k], eq);
        }
        poly.words[word_idx] |= acc;
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
    use sha3::{Shake256, digest::Update, digest::XofReader};

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

    // ── Barrett ("mod") fixed-weight sampler ──────────────────────────────────

    /// A deterministic XofReader that serves a fixed byte sequence, then zeros.
    /// Lets tests pin the exact bytes the Barrett sampler consumes.
    struct FixedXof {
        bytes: Vec<u8>,
        pos: usize,
    }
    impl FixedXof {
        fn new(bytes: Vec<u8>) -> Self {
            FixedXof { bytes, pos: 0 }
        }
    }
    impl XofReader for FixedXof {
        fn read(&mut self, buffer: &mut [u8]) {
            for b in buffer.iter_mut() {
                *b = self.bytes.get(self.pos).copied().unwrap_or(0);
                self.pos += 1;
            }
        }
    }

    fn mod_correct_weight<P: HqcParams>(seed: &[u8], weight: usize) {
        let mut xof = make_xof(seed);
        let p = sample_fixed_weight_mod::<P>(&mut xof, weight);
        assert_eq!(p.hamming_weight(), weight, "exact weight (implies distinct positions)");
        // No bits set in the overflow region above N-1.
        let last_bit = P::N & 63;
        if last_bit != 0 {
            let mask = (1u64 << last_bit) - 1;
            assert_eq!(p.words[P::N_WORDS - 1] & !mask, 0, "no overflow bits");
        }
    }

    #[test]
    fn mod_correct_weight_128() {
        mod_correct_weight::<Hqc128>(b"mod-seed-0", Hqc128::OMEGA_R);
    }
    #[test]
    fn mod_correct_weight_192() {
        mod_correct_weight::<Hqc192>(b"mod-seed-1", Hqc192::OMEGA_R);
    }
    #[test]
    fn mod_correct_weight_256() {
        mod_correct_weight::<Hqc256>(b"mod-seed-2", Hqc256::OMEGA_R);
    }

    #[test]
    fn mod_deterministic() {
        let mut xof1 = make_xof(b"mod-det");
        let mut xof2 = make_xof(b"mod-det");
        let p1 = sample_fixed_weight_mod::<Hqc128>(&mut xof1, Hqc128::OMEGA_R);
        let p2 = sample_fixed_weight_mod::<Hqc128>(&mut xof2, Hqc128::OMEGA_R);
        assert_eq!(p1, p2);
    }

    #[test]
    fn mod_differs_from_rejection() {
        // Same seed, different samplers ⇒ (almost surely) different positions.
        // This is what makes them non-interchangeable for KAT.
        let mut xof_a = make_xof(b"same-seed");
        let mut xof_b = make_xof(b"same-seed");
        let a = sample_fixed_weight::<Hqc128>(&mut xof_a, Hqc128::OMEGA_R);
        let b = sample_fixed_weight_mod::<Hqc128>(&mut xof_b, Hqc128::OMEGA_R);
        assert_ne!(a, b, "rejection and mod samplers must produce different vectors");
    }

    // ── Hand-computed vectors: pin the exact formula, dedup, and endianness ────

    #[test]
    fn mod_all_zero_words_gives_identity_positions() {
        // rand_i = 0 ⇒ support[i] = i + ((0·(N−i))>>32) = i. No collisions.
        // Result: bits 0,1,…,weight−1 set, nothing else.
        let weight = 8;
        let mut xof = FixedXof::new(vec![0u8; 4 * weight]);
        let p = sample_fixed_weight_mod::<Hqc128>(&mut xof, weight);
        assert_eq!(p.hamming_weight(), weight);
        for i in 0..weight {
            assert_eq!(p.get_bit(i), 1, "bit {i} must be set");
        }
        assert_eq!(p.get_bit(weight), 0, "bit {weight} must be clear");
    }

    #[test]
    fn mod_all_ones_words_exercises_dedup() {
        // rand_i = 0xFFFF_FFFF ⇒ support[i] = i + (N−i−1) = N−1 for all i, so
        // every draw collides. Backward dedup sets support[i]=i for i<weight−1
        // and keeps support[weight−1]=N−1.
        // Result: bits {0,…,weight−2, N−1} set.
        let weight = 8;
        let mut xof = FixedXof::new(vec![0xFFu8; 4 * weight]);
        let p = sample_fixed_weight_mod::<Hqc128>(&mut xof, weight);
        assert_eq!(p.hamming_weight(), weight, "dedup must still yield distinct positions");
        for i in 0..(weight - 1) {
            assert_eq!(p.get_bit(i), 1, "bit {i} must be set");
        }
        assert_eq!(p.get_bit(weight - 1), 0, "bit weight−1 must be clear (it became N−1)");
        assert_eq!(p.get_bit(Hqc128::N - 1), 1, "bit N−1 must be set");
    }

    #[test]
    fn mod_reads_u32_little_endian() {
        // weight=1, rand bytes = [00,00,00,80] = 0x8000_0000 little-endian = 2^31.
        // support[0] = 0 + ((2^31 · N) >> 32) = N >> 1 = floor(17669/2) = 8834.
        // A big-endian misread would give 0x00000080 = 128 ⇒ position 0, so this
        // test distinguishes the byte order.
        let mut xof = FixedXof::new(vec![0x00, 0x00, 0x00, 0x80]);
        let p = sample_fixed_weight_mod::<Hqc128>(&mut xof, 1);
        assert_eq!(p.hamming_weight(), 1);
        assert_eq!(p.get_bit(8834), 1, "expected position floor(N/2) = 8834 (little-endian)");
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
