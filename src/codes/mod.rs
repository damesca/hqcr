// Top-level concatenated codec: C.Encode and C.Decode (RMRS, spec §3.4).
//
// Encoding order: RS (outer) then RM (inner).
//   1. RS.Encode the k-byte message into n1 GF(2^8) symbols.
//   2. RM.Encode each symbol into one duplicated RM(1,7) block of N2 bits.
//   3. Concatenate the n1 blocks → n1 * N2 bits, embedded in the low bits of a
//      Poly<P> of ring dimension N. The trailing ℓ = N - n1*N2 bits stay zero.
//
// Decoding order: RM (inner) then RS (outer).
//   1. Split the low n1*N2 bits into n1 blocks of N2 bits.
//   2. RM.Decode each block (duplicated FHT) → one GF(2^8) symbol.
//   3. RS.Decode the n1 symbols → the k-byte message (or None on failure).
//
// ── Bit layout ────────────────────────────────────────────────────────────────
//   Both the RM byte buffers and Poly<P> use identical little-endian bit
//   packing: bit i lives at byte[i/8] position i%8 (equivalently word[i/64] at
//   i%64). Block j occupies bit range [j*N2, (j+1)*N2), which is byte-aligned
//   because N2 = 128 * MULTIPLICITY is a multiple of 8. This lets us build the
//   codeword in a flat byte buffer and copy it word-for-word into the Poly.

pub mod reed_muller;
pub mod reed_solomon;

use crate::params::HqcParams;
use crate::poly::Poly;

/// Bytes occupied by the full concatenated codeword: n1 blocks of N2 bits each.
/// N2 is a multiple of 8, so this division is exact.
#[inline]
fn codeword_bytes<P: HqcParams>() -> usize {
    P::N1 * P::N2 / 8
}

/// Bytes per RM block (one RS symbol after RM encoding): N2 / 8 = 16 * MULTIPLICITY.
#[inline]
fn block_bytes<P: HqcParams>() -> usize {
    P::N2 / 8
}

// ── C.Encode ──────────────────────────────────────────────────────────────────

/// Concatenated RMRS encode: `msg` (K bytes) → codeword embedded in `Poly<P>`.
///
/// The codeword occupies the low n1*N2 bits; the trailing ℓ = N - n1*N2 bits
/// remain zero. The result is added to `s·r2 + e` in PKE.Encrypt.
pub fn encode<P: HqcParams>(msg: &[u8]) -> Poly<P> {
    debug_assert_eq!(msg.len(), P::K, "message must be exactly K bytes");

    // Step 1: RS encode → n1 GF(2^8) symbols.
    let mut rs_cw = vec![0u8; P::N1];
    reed_solomon::rs_encode(msg, &mut rs_cw, P::DELTA);

    // Step 2 + 3: RM encode each symbol into its N2-bit block, concatenated.
    let mut buf = vec![0u8; codeword_bytes::<P>()];
    for (j, &symbol) in rs_cw.iter().enumerate() {
        reed_muller::rm_encode(symbol, P::MULTIPLICITY, &mut buf, j * P::N2);
    }

    // Pack the flat byte buffer into the low words of the Poly (little-endian).
    bytes_to_poly::<P>(&buf)
}

// ── C.Decode ──────────────────────────────────────────────────────────────────

/// Concatenated RMRS decode: `poly` (the noisy codeword `C.Encode(m) + err`)
/// → `Some(msg)` (K bytes) on success, or `None` if RS decoding fails.
pub fn decode<P: HqcParams>(poly: &Poly<P>) -> Option<Vec<u8>> {
    // Extract the low n1*N2 bits back into a flat byte buffer.
    let buf = poly_to_bytes::<P>(codeword_bytes::<P>(), poly);
    let bb = block_bytes::<P>();

    // Step 1 + 2: RM decode each N2-bit block → one GF(2^8) symbol.
    let mut rs_cw = vec![0u8; P::N1];
    for j in 0..P::N1 {
        let block = &buf[j * bb..(j + 1) * bb];
        rs_cw[j] = reed_muller::rm_decode(block, P::MULTIPLICITY);
    }

    // Step 3: RS decode the n1 symbols → message.
    reed_solomon::rs_decode(&rs_cw, P::K, P::DELTA)
}

// ── Poly ⟷ byte buffer (shared little-endian packing) ─────────────────────────

/// Pack a flat byte buffer into the low words of a fresh `Poly<P>`.
/// Byte k lands at word k/8, byte position k%8 (little-endian within the word).
fn bytes_to_poly<P: HqcParams>(buf: &[u8]) -> Poly<P> {
    debug_assert!(buf.len() <= P::N_WORDS * 8);
    let mut p = Poly::<P>::zero();
    for (k, &b) in buf.iter().enumerate() {
        p.words[k / 8] |= (b as u64) << (8 * (k % 8));
    }
    p
}

/// Extract the low `len` bytes from a `Poly<P>` (inverse of `bytes_to_poly`).
fn poly_to_bytes<P: HqcParams>(len: usize, p: &Poly<P>) -> Vec<u8> {
    debug_assert!(len <= P::N_WORDS * 8);
    let mut buf = vec![0u8; len];
    for (k, b) in buf.iter_mut().enumerate() {
        *b = (p.words[k / 8] >> (8 * (k % 8))) as u8;
    }
    buf
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::{Hqc128, Hqc192, Hqc256};

    /// Deterministic K-byte test message.
    fn test_msg<P: HqcParams>() -> Vec<u8> {
        (0..P::K)
            .map(|i| (i.wrapping_mul(37).wrapping_add(11)) as u8)
            .collect()
    }

    // ── Roundtrip with no errors ──────────────────────────────────────────────

    fn roundtrip_no_errors<P: HqcParams>() {
        let msg = test_msg::<P>();
        let cw = encode::<P>(&msg);
        let recovered = decode::<P>(&cw).expect("decode failed on clean codeword");
        assert_eq!(recovered, msg);
    }

    #[test]
    fn roundtrip_no_errors_128() {
        roundtrip_no_errors::<Hqc128>();
    }
    #[test]
    fn roundtrip_no_errors_192() {
        roundtrip_no_errors::<Hqc192>();
    }
    #[test]
    fn roundtrip_no_errors_256() {
        roundtrip_no_errors::<Hqc256>();
    }

    // ── Codeword occupies exactly the low n1*N2 bits ──────────────────────────

    #[test]
    fn codeword_fits_in_ring_with_zero_tail() {
        // The trailing ℓ = N - n1*N2 bits must be zero.
        let msg = test_msg::<Hqc128>();
        let cw = encode::<Hqc128>(&msg);
        for i in (Hqc128::N1 * Hqc128::N2)..Hqc128::N {
            assert_eq!(cw.get_bit(i), 0, "trailing bit {i} must be zero");
        }
    }

    // ── A few bit errors within RM capacity per block → recovers ──────────────

    fn flip_bit<P: HqcParams>(p: &mut Poly<P>, pos: usize) {
        if p.get_bit(pos) == 1 {
            p.clear_bit(pos);
        } else {
            p.set_bit(pos);
        }
    }

    fn roundtrip_small_bit_errors<P: HqcParams>() {
        // Flip 5 bits inside each of the first few blocks — far below RM's
        // per-block capacity (≥95 for ×3, ≥159 for ×5), so every symbol survives.
        let msg = test_msg::<P>();
        let mut cw = encode::<P>(&msg);
        for j in 0..P::N1.min(4) {
            for b in 0..5 {
                flip_bit::<P>(&mut cw, j * P::N2 + b * 7);
            }
        }
        let recovered = decode::<P>(&cw).expect("decode failed within RM capacity");
        assert_eq!(recovered, msg);
    }

    #[test]
    fn roundtrip_small_bit_errors_128() {
        roundtrip_small_bit_errors::<Hqc128>();
    }
    #[test]
    fn roundtrip_small_bit_errors_192() {
        roundtrip_small_bit_errors::<Hqc192>();
    }
    #[test]
    fn roundtrip_small_bit_errors_256() {
        roundtrip_small_bit_errors::<Hqc256>();
    }

    // ── δ fully-corrupted blocks → RM yields wrong symbols, RS fixes them ──────

    fn roundtrip_delta_symbol_errors<P: HqcParams>() {
        // Complement every bit of the first δ blocks. Flipping all bits of an
        // RM(1,7) codeword flips its all-ones bit, so RM decodes a *different*
        // symbol (MSB toggled). That is ≤ δ symbol errors, which RS corrects.
        let msg = test_msg::<P>();
        let mut cw = encode::<P>(&msg);
        for j in 0..P::DELTA {
            for b in 0..P::N2 {
                flip_bit::<P>(&mut cw, j * P::N2 + b);
            }
        }
        let recovered = decode::<P>(&cw).expect("RS should correct δ symbol errors");
        assert_eq!(recovered, msg);
    }

    #[test]
    fn roundtrip_delta_symbol_errors_128() {
        roundtrip_delta_symbol_errors::<Hqc128>();
    }
    #[test]
    fn roundtrip_delta_symbol_errors_192() {
        roundtrip_delta_symbol_errors::<Hqc192>();
    }
    #[test]
    fn roundtrip_delta_symbol_errors_256() {
        roundtrip_delta_symbol_errors::<Hqc256>();
    }
}
