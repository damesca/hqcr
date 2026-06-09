// HQC-PKE: the IND-CPA public-key encryption scheme underlying the KEM.
//
// Three operations (spec §4 / hqc_notes.md §4):
//
//   Keygen(seed_pke) -> (ek, dk)
//       (seed_dk, seed_ek) = I(seed_pke)           // SHA3-512 split, see note
//       (y, x) = sample_fixed_weight(XOF(seed_dk)) // sparse, weight ω
//       h      = sample_uniform(XOF(seed_ek))      // dense, public
//       s      = x + h·y
//       ek = (seed_ek, s)    dk = seed_dk
//
//   Encrypt(ek, m, θ) -> (u, v)
//       h           = sample_uniform(XOF(seed_ek)) // recomputed from ek
//       (r2, e, r1) = sample_fixed_weight(XOF(θ))  // sparse, weight ωr/ωe
//       u = r1 + h·r2
//       v = C.Encode(m) + s·r2 + e                 // truncated to n1·n2 bits
//
//   Decrypt(dk, (u, v)) -> Option<m>
//       y   = sample_fixed_weight(XOF(seed_dk))    // re-derived secret
//       tmp = v + u·y   = C.Encode(m) + (x·r2 + r1·y + e)
//       m   = C.Decode(tmp)                        // None on decode failure
//
// Sampling order matters and matches the reference exactly: y before x from
// seed_dk; r2, e, r1 (in that order) from θ. A KAT mismatch on ek/ct bytes
// would point here first.
//
// ── Architecture note: where the seed split lives ─────────────────────────────
// This follows the authoritative 2025 reference (Saarinen's hqc-py), NOT the
// simplified "Keygen(seed_dk, seed_ek)" sketch in CLAUDE.md. The KEM (step 12)
// expands seed_KEM via the XOF into (seed_pke, σ). One level down, here,
// PKE.Keygen expands seed_pke into (seed_dk, seed_ek) via the "I" function =
// SHA3-512(seed_pke ‖ 0x02), provided by `hash::i_pke_seed`. The XOF calls here
// go through `hash::xof`, which appends the 0x01 domain separator. Domain
// separators are fully wired (see hash.rs); the only remaining KAT gap is the
// sampler below.
//
// ── Sampler split (spec 2025) ─────────────────────────────────────────────────
// The two fixed-weight roles use different samplers, matching the reference:
//   • x, y (long-term secret) ← sample_fixed_weight       (rejection sampling)
//   • r2, e, r1 (ephemeral)   ← sample_fixed_weight_mod    (Barrett "mod" sampler)
// See poly/sampling.rs for why, and for the exact Barrett procedure.

use zeroize::{Zeroize, ZeroizeOnDrop};

use crate::codes;
use crate::hash;
use crate::params::{HqcParams, SEED_BYTES};
use crate::parsing;
use crate::poly::mul::{mul_dense_ct, mul_sparse_dense};
use crate::poly::sampling::{sample_fixed_weight, sample_fixed_weight_mod, sample_uniform};
use crate::poly::Poly;

// ── Key types ───────────────────────────────────────────────────────────────

/// HQC-PKE encryption key: the public seed for `h` plus the syndrome `s`.
/// Serializes (via `parsing`) to `seed_ek ‖ s`, identical to ekKEM.
pub struct EncryptionKey<P: HqcParams> {
    pub seed_ek: [u8; SEED_BYTES],
    pub s: Poly<P>,
}

// Manual Clone: the derived impl would add a spurious `where P: Clone` bound.
// `P` is a zero-sized marker, so cloning only copies `seed_ek` and `s`.
impl<P: HqcParams> Clone for EncryptionKey<P> {
    fn clone(&self) -> Self {
        Self {
            seed_ek: self.seed_ek,
            s: self.s.clone(),
        }
    }
}

/// HQC-PKE decryption key: just the 32-byte secret seed from which `y` is
/// re-derived on demand. Zeroized on drop.
#[derive(Zeroize, ZeroizeOnDrop)]
pub struct DecryptionKey {
    pub(crate) seed_dk: [u8; SEED_BYTES],
}

impl<P: HqcParams> EncryptionKey<P> {
    /// Serialize to the ekKEM wire format (`seed_ek ‖ s`), length `P::PK_BYTES`.
    pub fn to_bytes(&self) -> Vec<u8> {
        parsing::pack_public_key::<P>(&self.seed_ek, &self.s)
    }

    /// Parse from the ekKEM wire format. `None` on wrong length.
    pub fn from_bytes(bytes: &[u8]) -> Option<Self> {
        let (seed_ek, s) = parsing::unpack_public_key::<P>(bytes)?;
        Some(Self { seed_ek, s })
    }
}

// ── Keygen ──────────────────────────────────────────────────────────────────

/// HQC-PKE.Keygen. Deterministically derives the keypair from `seed_pke`.
pub fn keygen<P: HqcParams>(seed_pke: &[u8; SEED_BYTES]) -> (EncryptionKey<P>, DecryptionKey) {
    // I(seed_pke): SHA3-512 split into the dk-seed (secret) and ek-seed (public).
    let (seed_dk, seed_ek) = hash::i_pke_seed(seed_pke);

    // (y, x) ← sample_fixed_weight ×2 from XOF(seed_dk); y first, then x.
    let mut xof_dk = hash::xof(&seed_dk[..]);
    let y = sample_fixed_weight::<P>(&mut xof_dk, P::OMEGA);
    let x = sample_fixed_weight::<P>(&mut xof_dk, P::OMEGA);

    // h ← uniform from XOF(seed_ek).
    let mut xof_ek = hash::xof(&seed_ek[..]);
    let h = sample_uniform::<P>(&mut xof_ek);

    // s = x + h·y (y sparse ⇒ Mode-A multiply; positions are public here).
    let hy = mul_sparse_dense::<P>(&y, &h);
    let s = x.add(&hy);

    (
        EncryptionKey { seed_ek, s },
        DecryptionKey { seed_dk: *seed_dk },
    )
}

// ── Encrypt ─────────────────────────────────────────────────────────────────

/// HQC-PKE.Encrypt. Returns the ciphertext pair `(u, v)`.
///
/// `v` is truncated to its low `n1·n2` bits (the codeword region); the trailing
/// ring bits are zeroed so the Poly matches what is serialized on the wire.
pub fn encrypt<P: HqcParams>(ek: &EncryptionKey<P>, m: &[u8], theta: &[u8]) -> (Poly<P>, Poly<P>) {
    debug_assert_eq!(m.len(), P::K, "message must be exactly K bytes");

    // Recompute h from the public seed (the reference does not store h in ek).
    let mut xof_ek = hash::xof(&ek.seed_ek[..]);
    let h = sample_uniform::<P>(&mut xof_ek);

    // (r2, e, r1) ← sample_fixed_weight_mod ×3 from XOF(θ), in this exact order.
    // ωe = ωr, so all three use P::OMEGA_R. Ephemeral vectors use the Barrett
    // ("mod") sampler, not the rejection sampler used for the secret key.
    let mut xof_th = hash::xof(theta);
    let r2 = sample_fixed_weight_mod::<P>(&mut xof_th, P::OMEGA_R);
    let e = sample_fixed_weight_mod::<P>(&mut xof_th, P::OMEGA_R);
    let r1 = sample_fixed_weight_mod::<P>(&mut xof_th, P::OMEGA_R);

    // u = r1 + h·r2 (r2 sparse ⇒ Mode A).
    let hr2 = mul_sparse_dense::<P>(&r2, &h);
    let u = r1.add(&hr2);

    // v = C.Encode(m) + s·r2 + e, then truncate to the codeword bit-length.
    let cm = codes::encode::<P>(m);
    let sr2 = mul_sparse_dense::<P>(&r2, &ek.s);
    let mut v = cm.add(&sr2);
    v.add_assign(&e);
    truncate_to_bits::<P>(&mut v, P::N1 * P::N2);

    (u, v)
}

// ── Decrypt ─────────────────────────────────────────────────────────────────

/// HQC-PKE.Decrypt. Returns `Some(m)` (K bytes) on success, `None` if the inner
/// code fails to decode.
pub fn decrypt<P: HqcParams>(dk: &DecryptionKey, u: &Poly<P>, v: &Poly<P>) -> Option<Vec<u8>> {
    // Re-derive the secret y from seed_dk.
    let mut xof_dk = hash::xof(&dk.seed_dk[..]);
    let y = sample_fixed_weight::<P>(&mut xof_dk, P::OMEGA);

    // tmp = v + u·y. y is secret ⇒ Mode-B constant-time multiply on y.
    let uy = mul_dense_ct::<P>(u, &y);
    let tmp = v.add(&uy);

    codes::decode::<P>(&tmp)
}

// ── Helpers ───────────────────────────────────────────────────────────────────

/// Zero every bit at index `>= nbits` in `p` (clears the ring tail above the
/// codeword region). `nbits` is `n1·n2`, a multiple of 64 for all parameter
/// sets, but the partial-word path is kept for generality.
fn truncate_to_bits<P: HqcParams>(p: &mut Poly<P>, nbits: usize) {
    let full_words = nbits / 64;
    let rem = nbits % 64;
    if rem != 0 {
        p.words[full_words] &= (1u64 << rem) - 1;
        for w in &mut p.words[full_words + 1..P::N_WORDS] {
            *w = 0;
        }
    } else {
        for w in &mut p.words[full_words..P::N_WORDS] {
            *w = 0;
        }
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::{Hqc128, Hqc192, Hqc256};

    fn test_msg<P: HqcParams>() -> Vec<u8> {
        (0..P::K)
            .map(|i| (i.wrapping_mul(7).wrapping_add(1)) as u8)
            .collect()
    }

    // ── The core correctness oracle: Decrypt(Encrypt(m)) == m ─────────────────

    fn pke_roundtrip<P: HqcParams>() {
        let seed_pke = [0x42u8; SEED_BYTES];
        let (ek, dk) = keygen::<P>(&seed_pke);

        let m = test_msg::<P>();
        let theta = [0x17u8; 32];
        let (u, v) = encrypt::<P>(&ek, &m, &theta);

        let recovered = decrypt::<P>(&dk, &u, &v).expect("decrypt returned None");
        assert_eq!(recovered, m, "round-trip mismatch");
    }

    #[test]
    fn pke_roundtrip_128() {
        pke_roundtrip::<Hqc128>();
    }
    #[test]
    fn pke_roundtrip_192() {
        pke_roundtrip::<Hqc192>();
    }
    #[test]
    fn pke_roundtrip_256() {
        pke_roundtrip::<Hqc256>();
    }

    // ── Round-trip over several (seed, message, θ) triples ────────────────────

    #[test]
    fn pke_roundtrip_many_128() {
        for t in 0u8..8 {
            let seed_pke = [t.wrapping_mul(31).wrapping_add(3); SEED_BYTES];
            let (ek, dk) = keygen::<Hqc128>(&seed_pke);
            let m: Vec<u8> = (0..Hqc128::K)
                .map(|i| (i as u8) ^ t.wrapping_mul(5))
                .collect();
            let theta = [t ^ 0xA5; 32];
            let (u, v) = encrypt::<Hqc128>(&ek, &m, &theta);
            let got = decrypt::<Hqc128>(&dk, &u, &v).expect("decrypt None");
            assert_eq!(got, m, "mismatch at t={t}");
        }
    }

    // ── Determinism ───────────────────────────────────────────────────────────

    #[test]
    fn keygen_is_deterministic() {
        let seed = [0x5Au8; SEED_BYTES];
        let (ek1, dk1) = keygen::<Hqc128>(&seed);
        let (ek2, dk2) = keygen::<Hqc128>(&seed);
        assert_eq!(ek1.seed_ek, ek2.seed_ek);
        assert_eq!(ek1.s, ek2.s);
        assert_eq!(dk1.seed_dk, dk2.seed_dk);
    }

    #[test]
    fn encrypt_is_deterministic() {
        let (ek, _dk) = keygen::<Hqc128>(&[1u8; SEED_BYTES]);
        let m = test_msg::<Hqc128>();
        let theta = [9u8; 32];
        let (u1, v1) = encrypt::<Hqc128>(&ek, &m, &theta);
        let (u2, v2) = encrypt::<Hqc128>(&ek, &m, &theta);
        assert_eq!(u1, u2);
        assert_eq!(v1, v2);
    }

    // ── Different θ ⇒ different ciphertext, still decrypts to m ────────────────

    #[test]
    fn different_theta_changes_ciphertext() {
        let (ek, dk) = keygen::<Hqc128>(&[2u8; SEED_BYTES]);
        let m = test_msg::<Hqc128>();
        let (u1, v1) = encrypt::<Hqc128>(&ek, &m, &[1u8; 32]);
        let (u2, v2) = encrypt::<Hqc128>(&ek, &m, &[2u8; 32]);
        assert!(u1 != u2 || v1 != v2, "ciphertext should depend on θ");
        assert_eq!(decrypt::<Hqc128>(&dk, &u1, &v1).unwrap(), m);
        assert_eq!(decrypt::<Hqc128>(&dk, &u2, &v2).unwrap(), m);
    }

    // ── v has a zero tail above the codeword region ───────────────────────────

    #[test]
    fn encrypt_v_tail_is_zero() {
        let (ek, _dk) = keygen::<Hqc128>(&[3u8; SEED_BYTES]);
        let (_, v) = encrypt::<Hqc128>(&ek, &test_msg::<Hqc128>(), &[4u8; 32]);
        for i in (Hqc128::N1 * Hqc128::N2)..Hqc128::N {
            assert_eq!(v.get_bit(i), 0, "v bit {i} above codeword must be zero");
        }
    }

    // ── EncryptionKey byte round-trip ─────────────────────────────────────────

    #[test]
    fn encryption_key_byte_roundtrip() {
        let (ek, _dk) = keygen::<Hqc192>(&[7u8; SEED_BYTES]);
        let bytes = ek.to_bytes();
        assert_eq!(bytes.len(), Hqc192::PK_BYTES);
        let back = EncryptionKey::<Hqc192>::from_bytes(&bytes).expect("valid length");
        assert_eq!(back.seed_ek, ek.seed_ek);
        assert_eq!(back.s, ek.s);
    }
}
