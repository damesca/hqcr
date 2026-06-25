// Byte-level serialization matching the spec wire format exactly.
// The KAT vectors validate every byte produced here, so the bit-packing
// convention below is a hard correctness boundary.
//
// ── Bit-packing convention ────────────────────────────────────────────────────
//   A Poly<P> stores bit i at words[i/64] position i%64 (little-endian words).
//   On the wire a polynomial is a little-endian bit string: bit i lands at
//   byte[i/8] position i%8. These two layouts coincide, so packing is a plain
//   little-endian word→byte copy. The only subtlety is the final partial byte
//   when the bit count is not a multiple of 8 (true for ring elements, since N
//   is an odd prime): its high bits are forced to zero on pack and masked away
//   on unpack so the resulting Poly has no stray bits at or above the bit count.
//
// ── Wire formats (sizes asserted against params.rs) ───────────────────────────
//   ekKEM (public key):  seed_ek (32 B) || s     (⌈N/8⌉ B)            = PK_BYTES
//   cKEM  (ciphertext):  u (⌈N/8⌉ B) || v (⌈N1·N2/8⌉ B) || salt (16 B) = CT_BYTES
//   dkKEM (secret key):  compressed form is just seed_KEM (32 B) — no packing
//                        needed, so no function is provided here.
//
//   `s` and `u` are full ring elements (N bits). `v` is the concatenated RMRS
//   codeword (N1·N2 bits, a multiple of 8 because N2 = 128·MULTIPLICITY); its
//   trailing ℓ = N - N1·N2 ring bits are never transmitted. NOTE: N2 is the
//   *duplicated* RM length (384 or 640), so v is N1·N2 bits — not N1·128. This
//   matches params.rs CT_BYTES and the spec's published |cKEM| sizes.

use crate::params::{HqcParams, SALT_BYTES, SEED_BYTES};
use crate::poly::Poly;

// ── Byte-length helpers ───────────────────────────────────────────────────────

/// Bytes needed for a full ring element (N bits): ⌈N/8⌉.
#[inline]
pub fn ring_bytes<P: HqcParams>() -> usize {
    P::N.div_ceil(8)
}

/// Bytes needed for the codeword component v (N1·N2 bits, exact): N1·N2/8.
#[inline]
pub fn v_bytes<P: HqcParams>() -> usize {
    (P::N1 * P::N2).div_ceil(8)
}

// ── Core little-endian pack / unpack over an explicit bit count ────────────────

/// Pack the low `n_bits` of `poly` into ⌈n_bits/8⌉ little-endian bytes.
/// High bits of the final partial byte are zeroed.
fn pack<P: HqcParams>(poly: &Poly<P>, n_bits: usize) -> Vec<u8> {
    let n_bytes = n_bits.div_ceil(8);
    let mut out = vec![0u8; n_bytes];
    for (k, b) in out.iter_mut().enumerate() {
        *b = (poly.words[k / 8] >> (8 * (k % 8))) as u8;
    }
    // Zero the bits at/above n_bits in the final byte (no-op when byte-aligned).
    let rem = n_bits % 8;
    if rem != 0 {
        out[n_bytes - 1] &= (1u8 << rem) - 1;
    }
    out
}

/// Unpack ⌈n_bits/8⌉ little-endian bytes into a fresh `Poly<P>`.
/// Bits at/above `n_bits` are masked to zero so the Poly stays canonical.
fn unpack<P: HqcParams>(bytes: &[u8], n_bits: usize) -> Poly<P> {
    let n_bytes = n_bits.div_ceil(8);
    debug_assert_eq!(bytes.len(), n_bytes);
    let rem = n_bits % 8;

    let mut p = Poly::<P>::zero();
    for (k, &raw) in bytes.iter().enumerate() {
        // Strip stray high bits from the final partial byte before placing it.
        let byte = if rem != 0 && k == n_bytes - 1 {
            raw & ((1u8 << rem) - 1)
        } else {
            raw
        };
        p.words[k / 8] |= (byte as u64) << (8 * (k % 8));
    }
    p
}

// ── Ring element (s, u) ⟷ bytes ───────────────────────────────────────────────

/// Serialize a full ring element (N bits) to ⌈N/8⌉ bytes.
pub fn ring_to_bytes<P: HqcParams>(poly: &Poly<P>) -> Vec<u8> {
    pack(poly, P::N)
}

/// Deserialize ⌈N/8⌉ bytes into a ring element. Caller must pass exactly that many.
pub fn ring_from_bytes<P: HqcParams>(bytes: &[u8]) -> Poly<P> {
    unpack(bytes, P::N)
}

// ── Codeword component v ⟷ bytes ──────────────────────────────────────────────

/// Serialize the codeword component v (N1·N2 bits) to N1·N2/8 bytes.
/// Only the meaningful codeword bits are emitted; the ring tail is dropped.
pub fn v_to_bytes<P: HqcParams>(poly: &Poly<P>) -> Vec<u8> {
    pack(poly, P::N1 * P::N2)
}

/// Deserialize N1·N2/8 bytes into a ring element with the codeword in the low
/// bits and the trailing ℓ = N - N1·N2 bits zero.
pub fn v_from_bytes<P: HqcParams>(bytes: &[u8]) -> Poly<P> {
    unpack(bytes, P::N1 * P::N2)
}

// ── Public key: seed_ek (32 B) || s ───────────────────────────────────────────

/// Pack ekKEM = seed_ek || s. Output length is `P::PK_BYTES`.
pub fn pack_public_key<P: HqcParams>(seed_ek: &[u8; SEED_BYTES], s: &Poly<P>) -> Vec<u8> {
    let mut out = Vec::with_capacity(P::PK_BYTES);
    out.extend_from_slice(seed_ek);
    out.extend_from_slice(&ring_to_bytes(s));
    debug_assert_eq!(out.len(), P::PK_BYTES);
    out
}

/// Unpack ekKEM into (seed_ek, s). Returns `None` if the length is wrong.
pub fn unpack_public_key<P: HqcParams>(bytes: &[u8]) -> Option<([u8; SEED_BYTES], Poly<P>)> {
    if bytes.len() != P::PK_BYTES {
        return None;
    }
    let mut seed = [0u8; SEED_BYTES];
    seed.copy_from_slice(&bytes[..SEED_BYTES]);
    let s = ring_from_bytes::<P>(&bytes[SEED_BYTES..]);
    Some((seed, s))
}

// ── Ciphertext: u || v || salt ────────────────────────────────────────────────

/// Pack cKEM = u || v || salt. Output length is `P::CT_BYTES`.
pub fn pack_ciphertext<P: HqcParams>(u: &Poly<P>, v: &Poly<P>, salt: &[u8; SALT_BYTES]) -> Vec<u8> {
    let mut out = Vec::with_capacity(P::CT_BYTES);
    out.extend_from_slice(&ring_to_bytes(u));
    out.extend_from_slice(&v_to_bytes(v));
    out.extend_from_slice(salt);
    debug_assert_eq!(out.len(), P::CT_BYTES);
    out
}

/// Unpack cKEM into (u, v, salt). Returns `None` on any wrong length — this is
/// what lets Decaps reject a truncated ciphertext without panicking.
pub fn unpack_ciphertext<P: HqcParams>(
    bytes: &[u8],
) -> Option<(Poly<P>, Poly<P>, [u8; SALT_BYTES])> {
    if bytes.len() != P::CT_BYTES {
        return None;
    }
    let rb = ring_bytes::<P>();
    let vb = v_bytes::<P>();
    let u = ring_from_bytes::<P>(&bytes[..rb]);
    let v = v_from_bytes::<P>(&bytes[rb..rb + vb]);
    let mut salt = [0u8; SALT_BYTES];
    salt.copy_from_slice(&bytes[rb + vb..]);
    Some((u, v, salt))
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::{Hqc128, Hqc192, Hqc256};

    /// Deterministic poly with bits set at a fixed stride below `bound`.
    fn patterned_poly<P: HqcParams>(stride: usize, bound: usize) -> Poly<P> {
        let mut p = Poly::<P>::zero();
        let mut i = 1;
        while i < bound {
            p.set_bit(i);
            i += stride;
        }
        p
    }

    // ── Bit-order sanity: little-endian bit i at byte i/8 position i%8 ─────────

    #[test]
    fn bit_order_is_little_endian() {
        let mut p = Poly::<Hqc128>::zero();
        p.set_bit(0);
        p.set_bit(7);
        p.set_bit(8);
        p.set_bit(63);
        let bytes = ring_to_bytes(&p);
        assert_eq!(bytes[0], 0b1000_0001); // bits 0 and 7
        assert_eq!(bytes[1], 0b0000_0001); // bit 8
        assert_eq!(bytes[7], 0b1000_0000); // bit 63
    }

    // ── Output lengths match the spec constants ───────────────────────────────

    #[test]
    fn lengths_match_params() {
        assert_eq!(ring_bytes::<Hqc128>(), (Hqc128::N + 7) / 8);
        assert_eq!(
            ring_to_bytes(&Poly::<Hqc128>::zero()).len(),
            (Hqc128::N + 7) / 8
        );
        assert_eq!(
            v_to_bytes(&Poly::<Hqc128>::zero()).len(),
            Hqc128::N1 * Hqc128::N2 / 8
        );

        // Composite formats.
        let pk = pack_public_key::<Hqc128>(&[0u8; SEED_BYTES], &Poly::zero());
        assert_eq!(pk.len(), Hqc128::PK_BYTES);
        let ct = pack_ciphertext::<Hqc128>(&Poly::zero(), &Poly::zero(), &[0u8; SALT_BYTES]);
        assert_eq!(ct.len(), Hqc128::CT_BYTES);
    }

    fn check_lengths<P: HqcParams>(expected_pk: usize, expected_ct: usize) {
        assert_eq!(ring_bytes::<P>(), (P::N + 7) / 8);
        assert_eq!(ring_to_bytes(&Poly::<P>::zero()).len(), (P::N + 7) / 8);
        assert_eq!(v_to_bytes(&Poly::<P>::zero()).len(), P::N1 * P::N2 / 8);
        let pk = pack_public_key::<P>(&[0u8; SEED_BYTES], &Poly::zero());
        assert_eq!(pk.len(), P::PK_BYTES);
        assert_eq!(pk.len(), expected_pk);
        let ct = pack_ciphertext::<P>(&Poly::zero(), &Poly::zero(), &[0u8; SALT_BYTES]);
        assert_eq!(ct.len(), P::CT_BYTES);
        assert_eq!(ct.len(), expected_ct);
    }

    #[test]
    fn lengths_match_params_192() {
        check_lengths::<Hqc192>(4514, 8978);
    }
    #[test]
    fn lengths_match_params_256() {
        check_lengths::<Hqc256>(7237, 14421);
    }

    // ── Ring element round-trips for all three parameter sets ─────────────────

    fn ring_roundtrip<P: HqcParams>() {
        let p = patterned_poly::<P>(137, P::N);
        let bytes = ring_to_bytes(&p);
        assert_eq!(bytes.len(), ring_bytes::<P>());
        let back = ring_from_bytes::<P>(&bytes);
        assert_eq!(back, p);
    }

    #[test]
    fn ring_roundtrip_128() {
        ring_roundtrip::<Hqc128>();
    }
    #[test]
    fn ring_roundtrip_192() {
        ring_roundtrip::<Hqc192>();
    }
    #[test]
    fn ring_roundtrip_256() {
        ring_roundtrip::<Hqc256>();
    }

    // ── v round-trips (only the low N1*N2 bits survive) ───────────────────────

    fn v_roundtrip<P: HqcParams>() {
        let p = patterned_poly::<P>(91, P::N1 * P::N2);
        let bytes = v_to_bytes(&p);
        assert_eq!(bytes.len(), v_bytes::<P>());
        let back = v_from_bytes::<P>(&bytes);
        assert_eq!(back, p);
    }

    #[test]
    fn v_roundtrip_128() {
        v_roundtrip::<Hqc128>();
    }
    #[test]
    fn v_roundtrip_192() {
        v_roundtrip::<Hqc192>();
    }
    #[test]
    fn v_roundtrip_256() {
        v_roundtrip::<Hqc256>();
    }

    // ── v ring-tail bits (N1*N2 .. N-1) must not appear in v serialization ──────

    fn check_v_drops_ring_tail<P: HqcParams>() {
        // N1*N2 < N for all parameter sets; the gap is the ring tail that v drops.
        let tail_bit = P::N1 * P::N2;
        let base = patterned_poly::<P>(91, P::N1 * P::N2);
        let mut with_tail = base.clone();
        with_tail.set_bit(tail_bit);
        assert_eq!(
            v_to_bytes(&base),
            v_to_bytes(&with_tail),
            "ring-tail bit {tail_bit} must not appear in v serialization"
        );
    }

    #[test]
    fn v_drops_ring_tail_128() {
        check_v_drops_ring_tail::<Hqc128>();
    }
    #[test]
    fn v_drops_ring_tail_192() {
        check_v_drops_ring_tail::<Hqc192>();
    }
    #[test]
    fn v_drops_ring_tail_256() {
        check_v_drops_ring_tail::<Hqc256>();
    }

    // ── Partial final byte: high bits beyond N must be zero / masked ───────────

    #[test]
    fn final_byte_high_bits_are_zero_on_pack() {
        // Hqc128: N = 17669, N % 8 = 5, so the last byte has 3 high bits unused.
        let rem = Hqc128::N % 8;
        assert_ne!(rem, 0, "this test assumes N is not byte-aligned");
        let p = patterned_poly::<Hqc128>(3, Hqc128::N); // dense, exercises the tail
        let bytes = ring_to_bytes(&p);
        let last = *bytes.last().unwrap();
        assert_eq!(
            last & !((1u8 << rem) - 1),
            0,
            "high bits of final byte must be zero"
        );
    }

    fn check_final_byte_high_bits_zero<P: HqcParams>() {
        let rem = P::N % 8;
        assert_ne!(rem, 0, "this test assumes N is not byte-aligned");
        let p = patterned_poly::<P>(3, P::N);
        let bytes = ring_to_bytes(&p);
        let last = *bytes.last().unwrap();
        assert_eq!(last & !((1u8 << rem) - 1), 0, "high bits of final byte must be zero");
    }

    #[test]
    fn final_byte_high_bits_are_zero_192() {
        check_final_byte_high_bits_zero::<Hqc192>();
    }
    #[test]
    fn final_byte_high_bits_are_zero_256() {
        check_final_byte_high_bits_zero::<Hqc256>();
    }

    #[test]
    fn unpack_masks_stray_high_bits() {
        // Craft a ring buffer whose final byte has stray bits set above N%8.
        let mut bytes = ring_to_bytes(&patterned_poly::<Hqc128>(53, Hqc128::N));
        let rem = Hqc128::N % 8;
        *bytes.last_mut().unwrap() |= !((1u8 << rem) - 1); // set the unused high bits
        let p = ring_from_bytes::<Hqc128>(&bytes);
        // Re-packing must drop the stray bits, proving they never entered the Poly.
        let repacked = ring_to_bytes(&p);
        assert_eq!(repacked.last().unwrap() & !((1u8 << rem) - 1), 0);
    }

    fn check_unpack_masks_stray<P: HqcParams>() {
        let mut bytes = ring_to_bytes(&patterned_poly::<P>(53, P::N));
        let rem = P::N % 8;
        *bytes.last_mut().unwrap() |= !((1u8 << rem) - 1);
        let p = ring_from_bytes::<P>(&bytes);
        let repacked = ring_to_bytes(&p);
        assert_eq!(repacked.last().unwrap() & !((1u8 << rem) - 1), 0);
    }

    #[test]
    fn unpack_masks_stray_high_bits_192() {
        check_unpack_masks_stray::<Hqc192>();
    }
    #[test]
    fn unpack_masks_stray_high_bits_256() {
        check_unpack_masks_stray::<Hqc256>();
    }

    // ── Composite round-trips ─────────────────────────────────────────────────

    fn public_key_roundtrip<P: HqcParams>() {
        let seed: [u8; SEED_BYTES] =
            core::array::from_fn(|i| (i as u8).wrapping_mul(3).wrapping_add(1));
        let s = patterned_poly::<P>(101, P::N);
        let packed = pack_public_key::<P>(&seed, &s);
        let (seed_back, s_back) = unpack_public_key::<P>(&packed).expect("valid length");
        assert_eq!(seed_back, seed);
        assert_eq!(s_back, s);
    }

    #[test]
    fn public_key_roundtrip_128() {
        public_key_roundtrip::<Hqc128>();
    }
    #[test]
    fn public_key_roundtrip_192() {
        public_key_roundtrip::<Hqc192>();
    }
    #[test]
    fn public_key_roundtrip_256() {
        public_key_roundtrip::<Hqc256>();
    }

    fn ciphertext_roundtrip<P: HqcParams>() {
        let u = patterned_poly::<P>(71, P::N);
        let v = patterned_poly::<P>(83, P::N1 * P::N2);
        let salt: [u8; SALT_BYTES] = core::array::from_fn(|i| (i as u8) ^ 0x5A);
        let packed = pack_ciphertext::<P>(&u, &v, &salt);
        let (u_back, v_back, salt_back) = unpack_ciphertext::<P>(&packed).expect("valid length");
        assert_eq!(u_back, u);
        assert_eq!(v_back, v);
        assert_eq!(salt_back, salt);
    }

    #[test]
    fn ciphertext_roundtrip_128() {
        ciphertext_roundtrip::<Hqc128>();
    }
    #[test]
    fn ciphertext_roundtrip_192() {
        ciphertext_roundtrip::<Hqc192>();
    }
    #[test]
    fn ciphertext_roundtrip_256() {
        ciphertext_roundtrip::<Hqc256>();
    }

    // ── Wrong-length inputs are rejected, not panicked on ─────────────────────

    #[test]
    fn unpack_rejects_bad_lengths() {
        assert!(unpack_public_key::<Hqc128>(&[]).is_none());
        assert!(unpack_public_key::<Hqc128>(&vec![0u8; Hqc128::PK_BYTES - 1]).is_none());
        assert!(unpack_ciphertext::<Hqc128>(&[]).is_none());
        assert!(unpack_ciphertext::<Hqc128>(&vec![0u8; Hqc128::CT_BYTES + 1]).is_none());
    }
}
