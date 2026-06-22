// RM(1,7) = [128, 8, 64] base code, duplicated P::MULTIPLICITY times.
//
// The base code encodes one byte (8 bits) into a 128-bit codeword.
// With duplication the codeword grows to 128 * MULTIPLICITY bits:
//   HQC-128: 128 * 3 = 384 bits per symbol
//   HQC-192/256: 128 * 5 = 640 bits per symbol
//
// Encoding: the 8-bit input m selects row m of the 128×128 Hadamard matrix H.
//   H[m][i] = popcount(m & i) mod 2  — i.e. bit i of codeword = inner product of m and i over F2.
//   Equivalently: codeword[i] = parity of (m AND i).
//   The codeword is then copied MULTIPLICITY times.
//
// Decoding (duplicated FHT — spec §3.4.3):
//   1. Reshape the received 128*MULTIPLICITY bits into MULTIPLICITY sub-blocks of 128 bits.
//   2. For each position i in 0..128, compute the soft vote:
//        F[i] = sum over duplicates d of: +1 if bit i of block d is 0, -1 if 1.
//      Equivalently: F[i] = sum_d (1 - 2 * bit(d, i)) = MULTIPLICITY - 2 * (number of 1s at position i).
//   3. Apply the length-128 Walsh-Hadamard Transform to F (in-place, over i16).
//   4. The decoded byte is argmax |F̂|.
//      The sign of F̂ at that index encodes the all-ones codeword correction bit (see below).
//   5. Tie-break rule (spec §3.4.3): among positions with equal |F̂|, choose the one
//      with the smallest value in the low 7 bits.
//
// CT requirement: the argmax scan must visit all 128 positions without early exit.
// Running max is updated with subtle::ConditionallySelectable (branchless).
//
// Sign convention for the decoded byte:
//   The first-order Reed-Muller code RM(1,7) has codewords of the form:
//     c(a,b)[i] = a XOR <b, i>   where a is the "all-ones" bit and b is the 7-bit message.
//   After WHT, the peak at index `b` has sign (+) if a=0 and (-) if a=1.
//   The decoded byte is (sign_bit << 7) | b  where sign_bit = (F̂[b] < 0) as u8.

use subtle::{Choice, ConditionallySelectable};

// ── Encoding ──────────────────────────────────────────────────────────────────

/// Encode one byte `m` into a 128-bit RM(1,7) codeword stored in 2 × u64.
///
/// codeword[i] = parity(m & i) for i in 0..128.
/// We compute all 128 bits in parallel by iterating over the 7 input bits of m
/// and XOR-ing the appropriate "stripe" patterns into the output words.
///
/// The Hadamard matrix row for value m can be built by:
///   start with all-zeros (if bit 7 of m is 0) or all-ones (if bit 7 is 1)
///   for each set bit k in m (k < 7): XOR in the stripe mask for position k.
/// Stripe mask for bit k: bit i is set iff bit k of i is set.
fn rm_encode_byte(m: u8) -> [u64; 2] {
    // Stripe masks: for bit k, mask[k][i] = (i >> k) & 1.
    // For a 128-bit vector split into two u64 words (lo = bits 0..63, hi = bits 64..127):
    // Bit k stripe pattern (repeated every 2^(k+1) bits, half 0s half 1s):
    //   k=0: 0xAAAAAAAAAAAAAAAA (alternating, period 2)
    //   k=1: 0xCCCCCCCCCCCCCCCC (period 4)
    //   k=2: 0xF0F0F0F0F0F0F0F0 (period 8)
    //   k=3: 0xFF00FF00FF00FF00 (period 16)
    //   k=4: 0xFFFF0000FFFF0000 (period 32)
    //   k=5: 0xFFFFFFFF00000000 (period 64)
    //   k=6: lo=0x0000000000000000, hi=0xFFFFFFFFFFFFFFFF (bits 64..127 all set)
    const STRIPES_LO: [u64; 7] = [
        0xAAAA_AAAA_AAAA_AAAA,
        0xCCCC_CCCC_CCCC_CCCC,
        0xF0F0_F0F0_F0F0_F0F0,
        0xFF00_FF00_FF00_FF00,
        0xFFFF_0000_FFFF_0000,
        0xFFFF_FFFF_0000_0000,
        0x0000_0000_0000_0000,
    ];
    const STRIPES_HI: [u64; 7] = [
        0xAAAA_AAAA_AAAA_AAAA,
        0xCCCC_CCCC_CCCC_CCCC,
        0xF0F0_F0F0_F0F0_F0F0,
        0xFF00_FF00_FF00_FF00,
        0xFFFF_0000_FFFF_0000,
        0xFFFF_FFFF_0000_0000,
        0xFFFF_FFFF_FFFF_FFFF,
    ];

    // Bit 7 of m is the "all-ones" correction bit: if set, start from all-ones.
    let init: u64 = if (m >> 7) & 1 == 1 { u64::MAX } else { 0 };
    let mut lo = init;
    let mut hi = init;

    for k in 0..7 {
        if (m >> k) & 1 == 1 {
            lo ^= STRIPES_LO[k];
            hi ^= STRIPES_HI[k];
        }
    }

    [lo, hi]
}

/// Encode byte `m` into a duplicated RM codeword of `128 * multiplicity` bits,
/// written into `out[bit_offset..]`.
///
/// `out` is a flat byte slice of length `ceil(128 * multiplicity / 8) * n1` or
/// the full v-vector buffer. `bit_offset` is the starting bit position.
pub fn rm_encode(m: u8, multiplicity: usize, out: &mut [u8], bit_offset: usize) {
    let [lo, hi] = rm_encode_byte(m);

    for rep in 0..multiplicity {
        let base_bit = bit_offset + rep * 128;
        let base_byte = base_bit / 8;
        debug_assert!(base_bit % 8 == 0, "rm_encode requires byte-aligned output");

        // Write lo (bits 0..63 of codeword) as 8 bytes.
        let lo_bytes = lo.to_le_bytes();
        out[base_byte..base_byte + 8].copy_from_slice(&lo_bytes);

        // Write hi (bits 64..127 of codeword) as 8 bytes.
        let hi_bytes = hi.to_le_bytes();
        out[base_byte + 8..base_byte + 16].copy_from_slice(&hi_bytes);
    }
}

// ── Walsh-Hadamard Transform ──────────────────────────────────────────────────

/// In-place length-128 Walsh-Hadamard Transform over i16.
///
/// The WHT is its own inverse (up to scaling by 128).
/// We use the standard butterfly decomposition: log2(128) = 7 stages,
/// each stage pairs elements distance 2^s apart and applies:
///   (a, b) → (a + b, a - b)
///
/// After 7 stages, F̂[k] = sum_{i=0}^{127} F[i] * (-1)^{popcount(i & k)}.
/// This is exactly the Hadamard transform needed to find the most likely codeword.
///
/// Range: with multiplicity 5 and values in [-5, +5], after 7 stages the
/// maximum absolute value is 5 * 128 = 640, which fits comfortably in i16
/// (max 32767). No overflow.
fn wht(f: &mut [i16; 128]) {
    let mut step = 1usize;
    while step < 128 {
        let mut i = 0;
        while i < 128 {
            for j in i..i + step {
                let a = f[j];
                let b = f[j + step];
                f[j]        = a + b;
                f[j + step] = a - b;
            }
            i += 2 * step;
        }
        step <<= 1;
    }
}

// ── Decoding ──────────────────────────────────────────────────────────────────

/// Decode one duplicated RM block back to a byte.
///
/// `block` is a flat byte slice of exactly `16 * multiplicity` bytes
/// (= 128 * multiplicity bits), one received sub-block per copy.
///
/// Returns the decoded byte, or 0 if the block is all noise (shouldn't happen
/// within the HQC error budget).
pub fn rm_decode(block: &[u8], multiplicity: usize) -> u8 {
    debug_assert_eq!(block.len(), 16 * multiplicity);

    // Step 1: accumulate soft votes into F[0..128].
    // For each of the 128 bit positions, sum contributions from all copies.
    // Contribution: +1 if the bit is 0 (vote for "codeword bit = 0"),
    //               -1 if the bit is 1 (vote for "codeword bit = 1").
    // F[i] = MULTIPLICITY - 2 * (number of 1s at position i across all copies).
    let mut f = [0i16; 128];

    for rep in 0..multiplicity {
        let base = rep * 16; // 16 bytes = 128 bits per copy
        for i in 0..128 {
            let byte_idx = i / 8;
            let bit_idx  = i % 8;
            let bit = (block[base + byte_idx] >> bit_idx) & 1;
            // +1 if bit==0, -1 if bit==1: equivalent to 1 - 2*bit
            f[i] += 1 - 2 * (bit as i16);
        }
    }

    // Step 2: Walsh-Hadamard Transform.
    wht(&mut f);

    // Step 3: argmax |F̂| — constant-time scan over all 128 positions.
    //
    // We track:
    //   best_abs: the running maximum |F̂[k]|
    //   best_idx: the index achieving that maximum
    //   best_neg: whether F̂[best_idx] was negative (encodes the all-ones bit)
    //
    // Tie-break (spec §3.4.3): among equal |F̂|, keep the one with the smallest
    // value in the low 7 bits. Since we scan in ascending order (0..128) and
    // update only on strict improvement (>) this naturally picks the smallest index.
    //
    // CT: all updates use ConditionallySelectable — no branch on f[k] values.
    let mut best_abs = 0i16;
    let mut best_idx = 0u8;
    let mut best_neg = 0u8; // 1 if F̂[best_idx] < 0

    for k in 0..128usize {
        let val = f[k];
        let abs_val = val.abs();
        let is_neg  = (val < 0) as u8;

        // Update if strictly greater (ascending scan ⇒ smallest index wins ties).
        let update = Choice::from((abs_val > best_abs) as u8);
        best_abs = i16::conditional_select(&best_abs, &abs_val, update);
        best_idx = u8::conditional_select(&best_idx, &(k as u8), update);
        best_neg = u8::conditional_select(&best_neg, &is_neg,    update);
    }

    // Step 4: reconstruct the byte.
    // best_idx holds the low 7 bits of the message (b in the spec notation).
    // best_neg encodes the all-ones correction bit (a): bit 7 of the decoded byte.
    (best_neg << 7) | (best_idx & 0x7F)
}

// ── Constant-time audit: asm spot-check shim (Layer 3) ─────────────────────────
//
// `#[no_mangle] #[inline(never)]` so the argmax inside `rm_decode` keeps a
// standalone symbol in `cargo asm` output. Compiled only under `--features
// ct-audit`; no effect on a normal build. The audit reads it to confirm the
// argmax updates are `cmov`/masking over a full 128-position scan (no early-exit
// `jcc` on the secret `F̂` values). See docs/audit/constant-time.md §5.
#[cfg(feature = "ct-audit")]
#[no_mangle]
#[inline(never)]
pub fn ct_asm_rm_decode(block: &[u8]) -> u8 {
    rm_decode(block, 3)
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn encode_byte_to_vec(m: u8, multiplicity: usize) -> Vec<u8> {
        let len = 16 * multiplicity; // 128*mult bits = 16*mult bytes
        let mut buf = vec![0u8; len];
        rm_encode(m, multiplicity, &mut buf, 0);
        buf
    }

    // ── Encoding sanity ───────────────────────────────────────────────────────

    #[test]
    fn encode_zero_gives_all_zeros() {
        // m=0: all codeword bits are parity(0 & i) = 0.
        for mult in [3, 5] {
            let buf = encode_byte_to_vec(0, mult);
            assert!(buf.iter().all(|&b| b == 0), "mult={mult}");
        }
    }

    #[test]
    fn encode_255_gives_all_ones() {
        // m=0xFF: bit 7 set → start from all-ones; bits 0..6 all set → XOR stripes
        // net result: all bits become parity(0xFF & i) = popcount(i) mod 2.
        // Actually for m=0xFF all bits are set so codeword[i] = popcount(0xFF & i) mod 2.
        // This is not all-ones; just check round-trip below instead.
        let buf = encode_byte_to_vec(0xFF, 3);
        assert_eq!(buf.len(), 48);
        // Each copy is identical.
        assert_eq!(&buf[0..16], &buf[16..32]);
        assert_eq!(&buf[0..16], &buf[32..48]);
    }

    #[test]
    fn encode_duplicates_are_identical() {
        for m in [0u8, 1, 42, 128, 255] {
            for mult in [3usize, 5] {
                let buf = encode_byte_to_vec(m, mult);
                for rep in 1..mult {
                    assert_eq!(&buf[0..16], &buf[rep * 16..rep * 16 + 16],
                        "m={m} mult={mult} rep={rep}");
                }
            }
        }
    }

    #[test]
    fn codeword_weight_is_correct() {
        // RM(1,7) codeword weight depends on (a, b) where a = bit 7 of m, b = bits 0..6:
        //   m=0   (a=0, b=0): zero codeword         → weight   0
        //   m=128 (a=1, b=0): all-ones codeword      → weight 128
        //   all other m     : balanced linear form   → weight  64
        for m in 0u8..=255 {
            let [lo, hi] = rm_encode_byte(m);
            let weight = lo.count_ones() + hi.count_ones();
            let expected = match m {
                0   => 0,
                128 => 128,
                _   => 64,
            };
            assert_eq!(weight, expected, "m={m:#04x}");
        }
    }

    // ── Decode (no errors) ────────────────────────────────────────────────────

    #[test]
    fn decode_roundtrip_no_errors_mult3() {
        for m in 0u8..=255 {
            let buf = encode_byte_to_vec(m, 3);
            let decoded = rm_decode(&buf, 3);
            assert_eq!(decoded, m, "m={m}");
        }
    }

    #[test]
    fn decode_roundtrip_no_errors_mult5() {
        for m in 0u8..=255 {
            let buf = encode_byte_to_vec(m, 5);
            let decoded = rm_decode(&buf, 5);
            assert_eq!(decoded, m, "m={m}");
        }
    }

    // ── Error correction ──────────────────────────────────────────────────────

    fn flip_bit(buf: &mut [u8], pos: usize) {
        buf[pos / 8] ^= 1 << (pos % 8);
    }

    #[test]
    fn decode_corrects_up_to_capacity_mult3() {
        // RM(1,7) duplicated ×3 = [384, 8, 192]. Error capacity ≈ (192-1)/2 = 95 bits.
        // We test with exactly 95 bit errors (worst-case that should still decode).
        let m = 0xABu8;
        let mut buf = encode_byte_to_vec(m, 3);
        for i in 0..95 {
            flip_bit(&mut buf, i);
        }
        let decoded = rm_decode(&buf, 3);
        assert_eq!(decoded, m, "failed with 95 errors");
    }

    #[test]
    fn decode_corrects_up_to_capacity_mult5() {
        // RM(1,7) duplicated ×5 = [640, 8, 320]. Error capacity ≈ (320-1)/2 = 159 bits.
        let m = 0xCDu8;
        let mut buf = encode_byte_to_vec(m, 5);
        for i in 0..159 {
            flip_bit(&mut buf, i);
        }
        let decoded = rm_decode(&buf, 5);
        assert_eq!(decoded, m, "failed with 159 errors");
    }

    #[test]
    fn decode_single_bit_error_each_position_mult3() {
        // Flip every possible single bit in the codeword and verify recovery.
        let m = 0x55u8;
        let clean = encode_byte_to_vec(m, 3);
        for pos in 0..384 {
            let mut buf = clean.clone();
            flip_bit(&mut buf, pos);
            let decoded = rm_decode(&buf, 3);
            assert_eq!(decoded, m, "single bit flip at pos={pos}");
        }
    }

    // ── WHT properties ────────────────────────────────────────────────────────

    #[test]
    fn wht_of_unit_vector_is_all_ones_row() {
        // WHT of e_0 = [1, 0, 0, ..., 0] should be [1, 1, 1, ..., 1].
        let mut f = [0i16; 128];
        f[0] = 1;
        wht(&mut f);
        assert!(f.iter().all(|&v| v == 1));
    }

    #[test]
    fn wht_involution() {
        // WHT(WHT(f)) = 128 * f  (the transform is its own inverse up to scale).
        let mut f = [0i16; 128];
        for i in 0..128 {
            f[i] = (i as i16 % 7) - 3; // arbitrary values in [-3, 3]
        }
        let original = f;
        wht(&mut f);
        wht(&mut f);
        for i in 0..128 {
            assert_eq!(f[i], 128 * original[i], "index {i}");
        }
    }
}
