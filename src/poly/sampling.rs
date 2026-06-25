// Sampling functions for polynomials in R = F2[X]/(X^n - 1).
//
// Two fixed-weight samplers + one uniform, matching the 2025 reference
// (reference_impl/src/vector.c) exactly. The reference deliberately uses a
// DIFFERENT fixed-weight routine for the secret vs the ephemeral vectors, and
// hqc.c picks them per role — using the wrong one is a KAT mismatch:
//
//   sample_fixed_weight     — reference `vect_sample_fixed_weight1`
//                             (`vect_generate_random_support1`). Used for the
//                             SECRET x, y in keygen. Draws 24-bit BIG-ENDIAN
//                             values in 3·weight-byte batches, rejects any value
//                             ≥ ⌊2^24/N⌋·N (redraw), reduces survivors mod N
//                             (Barrett), and redraws on a duplicate position.
//   sample_fixed_weight_mod — reference `vect_sample_fixed_weight2`
//                             (`vect_generate_random_support2`). Used for the
//                             EPHEMERAL r1, r2, e in encrypt. Reads 4·weight
//                             bytes as little-endian u32 and sets
//                             support[i] = i + ((rand·(N−i))>>32), then a
//                             backward dedup; no rejection.
//   sample_uniform          — fills all N bits uniformly (public key h).
//
// History: an earlier revision used a single bespoke rejection sampler (u16) for
// x/y and the mod sampler for everything; reading the reference proved x/y need
// sampler #1 specifically. See KAT.md "Sampler fix".
//
// Constant-time contract:
//   sample_fixed_weight is REJECTION-based and therefore NOT constant time — its
//   running time and XOF consumption depend on the drawn values. This matches
//   the reference, which accepts variable timing for the once-per-keygen secret
//   sampling.
//   sample_fixed_weight_mod is branchless: the reduction, the duplicate
//   resolution, and the final bit-setting all run over fixed-length loops with
//   constant-time selects, leaking nothing about the sampled positions.

use sha3::digest::XofReader;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq};

use super::Poly;
use crate::params::HqcParams;

// ── Internal XOF helpers ──────────────────────────────────────────────────────

/// Read `out.len()` bytes from the XOF, then consume the trailing alignment
/// padding so the XOF advances by `ceil(len/8)·8` bytes — exactly matching the
/// reference `xof_get_bytes` (reference_impl/src/symmetric.c), which squeezes in
/// 8-byte (`sizeof(uint64_t)`) units and discards the unused tail of the final
/// unit. This inter-call discard is what keeps successive samples drawn from the
/// SAME XOF (e.g. x after y, or e/r1 after r2) byte-aligned with the reference.
fn xof_get_bytes(xof: &mut impl XofReader, out: &mut [u8]) {
    xof.read(out);
    let pad = (8 - out.len() % 8) % 8;
    if pad != 0 {
        let mut discard = [0u8; 8];
        xof.read(&mut discard[..pad]);
    }
}

/// Read exactly 8 bytes from the XOF and interpret them as a little-endian u64.
#[inline(always)]
fn read_u64(xof: &mut impl XofReader) -> u64 {
    let mut buf = [0u8; 8];
    xof.read(&mut buf);
    u64::from_le_bytes(buf)
}

// ── Fixed-weight sampler #1 (secret x, y) ─────────────────────────────────────

/// Sample a polynomial with exactly `weight` coefficients equal to 1, using the
/// reference `vect_sample_fixed_weight1` / `vect_generate_random_support1`
/// (reference_impl/src/vector.c). This is the SECRET-key sampler, for x and y.
///
/// Algorithm:
///   - Draw 24-bit BIG-ENDIAN values, refilling `xof` in `3·weight`-byte
///     batches (reading full batches reproduces the reference's byte
///     consumption — including the discard of a batch's unused tail — which is
///     what keeps a second call on the same XOF, x after y, byte-aligned).
///   - Reject any draw ≥ the rejection threshold `⌊2^24/N⌋·N` (redraw): this is
///     what makes the surviving values reduce *uniformly* mod N.
///   - Reduce the survivor mod N (the reference's constant-time `barrett_reduce`;
///     for x < 2^24 this equals `x % N`, which we use directly).
///   - Redraw on a duplicate position (do not advance).
///
/// NOT constant time: like the reference, the number of rejections / duplicate
/// redraws — and hence the running time and XOF byte consumption — depends on
/// the drawn values. The reference accepts this for the once-per-keygen secret
/// sampling. (Reproducing the reference's exact rejection behaviour is required
/// for KAT correctness, so this timing property is inherent, not a choice.)
pub fn sample_fixed_weight<P: HqcParams>(xof: &mut impl XofReader, weight: usize) -> Poly<P> {
    debug_assert!(
        weight <= 256,
        "weight {weight} exceeds internal buffer size"
    );
    debug_assert!(weight <= P::N, "weight {weight} > N={}", P::N);

    let n = P::N as u32;
    // Rejection threshold t = ⌊2^24 / N⌋ · N (UTILS_REJECTION_THRESHOLD).
    let threshold = ((1u32 << 24) / n) * n;

    // Random bytes are consumed in batches of `3·weight`, matching the
    // reference's `rand_bytes[3·weight]` refill granularity.
    let batch = 3 * weight;
    let mut buf = [0u8; 3 * 256];
    let mut j = batch; // == batch forces an initial refill on first use

    let mut positions = [0u32; 256];
    let mut filled: usize = 0;

    while filled < weight {
        // Draw a 24-bit big-endian candidate, rejecting values ≥ threshold.
        let pos = loop {
            if j == batch {
                // One reference `xof_get_bytes(3·weight)` call per batch: reads
                // `batch` bytes and discards the 8-byte-alignment tail.
                xof_get_bytes(xof, &mut buf[..batch]);
                j = 0;
            }
            let cand = ((buf[j] as u32) << 16) | ((buf[j + 1] as u32) << 8) | (buf[j + 2] as u32);
            j += 3;
            if cand < threshold {
                break cand % n; // == barrett_reduce(cand)
            }
        };

        positions[filled] = pos;

        // Redraw on duplicate: scan the accepted prefix; only advance if distinct.
        let mut dup = false;
        for &p in &positions[..filled] {
            if p == pos {
                dup = true;
            }
        }
        if !dup {
            filled += 1;
        }
    }

    // Build the polynomial: set bit at each accepted position.
    let mut poly = Poly::<P>::zero();
    for &p in &positions[..weight] {
        poly.set_bit(p as usize);
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
///      support[i] = i + ((rand_i · (N − i)) >> 32)
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
pub fn sample_fixed_weight_mod<P: HqcParams>(xof: &mut impl XofReader, weight: usize) -> Poly<P> {
    debug_assert!(
        weight <= 256,
        "weight {weight} exceeds internal buffer size"
    );
    debug_assert!(weight <= P::N, "weight {weight} > N={}", P::N);

    let mut support = [0u32; 256];

    // Draw all 4·weight random bytes in one reference `xof_get_bytes(4·weight)`
    // call (8-byte-aligned consumption with tail discard), then read them as
    // little-endian u32 — matching `vect_generate_random_support2`.
    let nbytes = 4 * weight;
    let mut buf = [0u8; 4 * 256];
    xof_get_bytes(xof, &mut buf[..nbytes]);

    // Step 1 — draw and Barrett-reduce: support[i] = i + ((rand_i·(N−i)) >> 32).
    for i in 0..weight {
        let rand =
            u32::from_le_bytes([buf[4 * i], buf[4 * i + 1], buf[4 * i + 2], buf[4 * i + 3]]) as u64;
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
    use sha3::{digest::Update, digest::XofReader, Shake256};

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

    fn check_fixed_weight_bits_in_range<P: HqcParams>(seed: &[u8]) {
        let mut xof = make_xof(seed);
        let p = sample_fixed_weight::<P>(&mut xof, P::OMEGA_R);
        let last_bit = P::N & 63;
        let mask = (1u64 << last_bit) - 1;
        assert_eq!(p.words[P::N_WORDS - 1] & !mask, 0, "no overflow bits");
    }

    #[test]
    fn fixed_weight_all_bits_in_range_128() {
        check_fixed_weight_bits_in_range::<Hqc128>(b"test-seed-range");
    }
    #[test]
    fn fixed_weight_all_bits_in_range_192() {
        check_fixed_weight_bits_in_range::<Hqc192>(b"test-seed-range-192");
    }
    #[test]
    fn fixed_weight_all_bits_in_range_256() {
        check_fixed_weight_bits_in_range::<Hqc256>(b"test-seed-range-256");
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
        assert_eq!(
            p.hamming_weight(),
            weight,
            "exact weight (implies distinct positions)"
        );
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
        assert_ne!(
            a, b,
            "rejection and mod samplers must produce different vectors"
        );
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
        assert_eq!(
            p.hamming_weight(),
            weight,
            "dedup must still yield distinct positions"
        );
        for i in 0..(weight - 1) {
            assert_eq!(p.get_bit(i), 1, "bit {i} must be set");
        }
        assert_eq!(
            p.get_bit(weight - 1),
            0,
            "bit weight−1 must be clear (it became N−1)"
        );
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
        assert_eq!(
            p.get_bit(8834),
            1,
            "expected position floor(N/2) = 8834 (little-endian)"
        );
    }

    // ── uniform sampler ───────────────────────────────────────────────────────

    fn check_uniform_no_overflow<P: HqcParams>(seed: &[u8]) {
        let mut xof = make_xof(seed);
        let p = sample_uniform::<P>(&mut xof);
        let last_bit = P::N & 63;
        let mask = (1u64 << last_bit) - 1;
        assert_eq!(p.words[P::N_WORDS - 1] & !mask, 0, "no overflow bits");
    }

    #[test]
    fn uniform_no_overflow_bits_128() {
        check_uniform_no_overflow::<Hqc128>(b"uniform-test");
    }
    #[test]
    fn uniform_no_overflow_bits_192() {
        check_uniform_no_overflow::<Hqc192>(b"uniform-192");
    }
    #[test]
    fn uniform_no_overflow_bits_256() {
        check_uniform_no_overflow::<Hqc256>(b"uniform-256");
    }

    #[test]
    fn uniform_deterministic() {
        let mut xof1 = make_xof(b"uniform-det");
        let mut xof2 = make_xof(b"uniform-det");
        let p1 = sample_uniform::<Hqc128>(&mut xof1);
        let p2 = sample_uniform::<Hqc128>(&mut xof2);
        assert_eq!(p1, p2);
    }
}
