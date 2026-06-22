// HQC-KEM: the IND-CCA2 key-encapsulation mechanism built on top of HQC-PKE
// via the salted Fujisaki-Okamoto transform with implicit rejection (SFO⊥_m,
// spec 2025 §5). Three operations:
//
//   Keygen(seed_KEM) -> (ek, dk)
//       (seed_pke, σ) = XOF(seed_KEM)              // SHAKE256 expansion
//       (ek, dk_PKE)  = PKE.Keygen(seed_pke)
//       ek = ekKEM (seed_ek ‖ s)   dk = (ek, dk_PKE, σ, seed_KEM)
//
//   Encaps(ek) -> (K, c)
//       m    ←$ random(k B)
//       salt ←$ random(16 B)
//       (K, θ) = G(H(ek), m, salt)
//       (u, v) = PKE.Encrypt(ek, m, θ)
//       c = u ‖ v ‖ salt
//
//   Decaps(dk, c) -> K
//       parse c = (u, v, salt)                     // None ⇒ implicit reject
//       m'      = PKE.Decrypt(dk_PKE, u, v)        // Option
//       (K', θ')= G(H(ek), m' or zeros, salt)
//       c'      = PKE.Encrypt(ek, m' or zeros, θ')
//       K̄       = J(H(ek), σ, c)
//       valid   = (m' ≠ ⊥) ∧ (c' == c)             // constant-time
//       return  ct_select(valid, K', K̄)
//
// ── Seed / σ derivation (Saarinen-style, matches pke.rs) ──────────────────────
// The compressed secret key is just the 32-byte seed_KEM, so σ MUST be
// derivable from it: we expand seed_KEM through the SHAKE256 XOF into
// (seed_pke ‖ σ), seed_pke first. PKE.Keygen then does its own SHA3-512 split
// of seed_pke into (seed_dk, seed_ek). See the architecture note in pke.rs.
//
// ── Constant-time discipline (spec §5, CLAUDE.md) ─────────────────────────────
// Decaps must not leak, via timing or the returned key, whether decryption
// succeeded. We therefore:
//   • always run G and the re-encryption, feeding an all-zero message on
//     decode failure (never branch the crypto on m' = ⊥);
//   • compare c' to c with `subtle::ConstantTimeEq`;
//   • fold the decode-success bit into a `subtle::Choice`;
//   • select K' vs K̄ with `subtle::ConditionallySelectable`.
// The only public-data branch is the ciphertext length check: a wrong-length c
// is rejected up front (its length carries no secret), returning J(H(ek), σ, c).

use rand_core::{CryptoRng, RngCore};
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq};
use zeroize::{Zeroize, Zeroizing};

use crate::hash::{self, SharedKey};
use crate::params::{HqcParams, SALT_BYTES, SEED_BYTES, SHARED_KEY_BYTES};
use crate::parsing;
use crate::pke::{self, DecryptionKey, EncryptionKey};

// ── Key types ───────────────────────────────────────────────────────────────

/// KEM encapsulation key (public). Identical to the PKE encryption key and to
/// the `ekKEM` wire format (`seed_ek ‖ s`).
pub type PublicKey<P> = EncryptionKey<P>;

/// KEM decapsulation key. The wire form is the compressed 32-byte `seed_KEM`
/// (`to_bytes`); everything below is re-derived from it (`from_bytes`).
///
/// Holds the secret σ and the PKE secret seed (both zeroized on drop) plus the
/// public encapsulation key, which Decaps needs to recompute `H(ek)` and to
/// re-encrypt.
pub struct DecapsulationKey<P: HqcParams> {
    seed_kem: Zeroizing<[u8; SEED_BYTES]>,
    sigma: Zeroizing<[u8; SEED_BYTES]>,
    dk_pke: DecryptionKey,
    ek: EncryptionKey<P>,
}

impl<P: HqcParams> DecapsulationKey<P> {
    /// The matching public encapsulation key.
    pub fn public_key(&self) -> &EncryptionKey<P> {
        &self.ek
    }

    /// Compressed secret key: the 32-byte `seed_KEM`. Zeroized on drop.
    pub fn to_bytes(&self) -> Zeroizing<[u8; SEED_BYTES]> {
        self.seed_kem.clone()
    }

    /// Reconstruct the full key from the compressed `seed_KEM`. `None` on wrong
    /// length. Equivalent to re-running `keygen_from_seed`.
    pub fn from_bytes(bytes: &[u8]) -> Option<Self> {
        if bytes.len() != SEED_BYTES {
            return None;
        }
        // Wipe this secret-seed copy on the way out — it is the compressed
        // secret key; `keygen_from_seed` only borrows it. (Zeroize audit G1.)
        let mut seed = Zeroizing::new([0u8; SEED_BYTES]);
        seed.copy_from_slice(bytes);
        let (_, dk) = keygen_from_seed::<P>(&seed);
        Some(dk)
    }

    /// **Debug / KAT only** (gated behind the `kat` feature). Serializes the
    /// decapsulation key in the HQC *reference* secret-key wire format used by
    /// the official NIST KAT `.rsp` files, rather than the compressed 32-byte
    /// `seed_KEM` returned by [`to_bytes`](Self::to_bytes).
    ///
    /// Matches the reference `crypto_kem_keypair` `dk_kem` layout (kem.c):
    ///
    /// ```text
    /// sk = ek_pke (PK_BYTES) ‖ dk_pke (32) ‖ σ (K) ‖ seed_kem (32)
    /// ```
    ///
    /// where `dk_pke = seed_dk` (the first half of the I-split, hqc.c
    /// `hqc_pke_keygen`), σ is the reference's `K`-byte randomness (this crate
    /// squeezes a 32-byte σ from the XOF, of which the first `K` bytes coincide
    /// with the reference), and the suffix `seed_kem` is the master KEM seed.
    /// Total suffix is `64 + K` bytes, matching the upstream `|sk|` of
    /// 2321 / 4602 / 7333.
    ///
    /// Exposed purely so the KAT harness can emit a reference-shaped `sk` for
    /// byte-for-byte comparison against pqc-hqc.org. NOT part of the production
    /// API and not constant-time-audited for this layout.
    #[cfg(feature = "kat")]
    pub fn expanded_secret_key_bytes(&self) -> Vec<u8> {
        let mut out = Vec::with_capacity(P::PK_BYTES + 2 * SEED_BYTES + P::K);
        out.extend_from_slice(&self.ek.to_bytes()); // ek_pke = seed_ek ‖ s
        out.extend_from_slice(&self.dk_pke.seed_dk); // dk_pke = seed_dk, 32 B
        out.extend_from_slice(&self.sigma[..P::K]); // σ, K B (reference length)
        out.extend_from_slice(&self.seed_kem[..]); // seed_kem, 32 B
        debug_assert_eq!(out.len(), P::PK_BYTES + 2 * SEED_BYTES + P::K);
        out
    }
}

// ── Keygen ──────────────────────────────────────────────────────────────────

/// HQC-KEM.Keygen with caller-supplied randomness. Draws a fresh 32-byte
/// `seed_KEM` from `rng` and expands it deterministically.
pub fn keygen<P: HqcParams, R: RngCore + CryptoRng>(
    rng: &mut R,
) -> (PublicKey<P>, DecapsulationKey<P>) {
    let mut seed_kem = Zeroizing::new([0u8; SEED_BYTES]);
    rng.fill_bytes(&mut seed_kem[..]);
    keygen_from_seed::<P>(&seed_kem)
}

/// Deterministic HQC-KEM.Keygen from a fixed `seed_KEM`. The compressed secret
/// key is exactly this seed; the rest is derived here.
pub fn keygen_from_seed<P: HqcParams>(
    seed_kem: &[u8; SEED_BYTES],
) -> (PublicKey<P>, DecapsulationKey<P>) {
    // Expand seed_KEM → (seed_pke ‖ σ) via SHAKE256 (seed_pke first).
    let mut xof = hash::xof(&seed_kem[..]);
    let mut seed_pke = Zeroizing::new([0u8; SEED_BYTES]);
    read_xof(&mut xof, &mut seed_pke[..]);
    let mut sigma = Zeroizing::new([0u8; SEED_BYTES]);
    read_xof(&mut xof, &mut sigma[..]);

    let (ek, dk_pke) = pke::keygen::<P>(&seed_pke);

    let dk = DecapsulationKey {
        seed_kem: Zeroizing::new(*seed_kem),
        sigma,
        dk_pke,
        ek: ek.clone(),
    };
    (ek, dk)
}

// ── Encaps ──────────────────────────────────────────────────────────────────

/// HQC-KEM.Encaps with caller-supplied randomness. Returns the shared key `K`
/// and the ciphertext `c` (length `P::CT_BYTES`).
pub fn encaps<P: HqcParams, R: RngCore + CryptoRng>(
    rng: &mut R,
    ek: &PublicKey<P>,
) -> (SharedKey, Vec<u8>) {
    let mut m = Zeroizing::new(vec![0u8; P::K]);
    rng.fill_bytes(&mut m[..]);
    let mut salt = [0u8; SALT_BYTES];
    rng.fill_bytes(&mut salt);
    encaps_deterministic::<P>(ek, &m, &salt)
}

/// Deterministic HQC-KEM.Encaps core: fixed message `m` (K bytes) and `salt`.
/// Exposed for testing and KAT reproduction.
pub fn encaps_deterministic<P: HqcParams>(
    ek: &PublicKey<P>,
    m: &[u8],
    salt: &[u8; SALT_BYTES],
) -> (SharedKey, Vec<u8>) {
    debug_assert_eq!(m.len(), P::K, "message must be exactly K bytes");

    let ek_bytes = ek.to_bytes();
    let ek_hash = hash::h_ek(&ek_bytes);

    // (K, θ) = G(H(ek), m, salt).
    let (k, theta) = hash::g(&ek_hash, m, salt);

    // (u, v) = PKE.Encrypt(ek, m, θ); ciphertext is u ‖ v ‖ salt.
    let (u, v) = pke::encrypt::<P>(ek, m, &theta[..]);
    let c = parsing::pack_ciphertext::<P>(&u, &v, salt);

    (k, c)
}

// ── Decaps ──────────────────────────────────────────────────────────────────

/// HQC-KEM.Decaps. Always returns a 32-byte key and never panics: a malformed,
/// truncated, or tampered ciphertext yields the implicit-rejection key derived
/// from σ, indistinguishable (by timing or value) from a valid decapsulation.
pub fn decaps<P: HqcParams>(dk: &DecapsulationKey<P>, c: &[u8]) -> SharedKey {
    let ek_bytes = dk.ek.to_bytes();
    let ek_hash = hash::h_ek(&ek_bytes);

    // Implicit-rejection key K̄ = J(H(ek), σ, c) — computed over the raw c.
    let mut k_bar = hash::j(&ek_hash, &dk.sigma, c);

    // A wrong-length ciphertext carries no secret; reject it up front with K̄.
    let (u, v, salt) = match parsing::unpack_ciphertext::<P>(c) {
        Some(parts) => parts,
        None => return k_bar,
    };

    // m' = PKE.Decrypt(dk_PKE, u, v). On failure, feed G/Encrypt all zeros so
    // the timing of both is independent of decode success.
    let m_prime = pke::decrypt::<P>(&dk.dk_pke, &u, &v);
    let decode_ok = Choice::from(m_prime.is_some() as u8);
    let m_bytes: Zeroizing<Vec<u8>> =
        Zeroizing::new(m_prime.unwrap_or_else(|| vec![0u8; P::K]));

    // Re-derive (K', θ') and re-encrypt under the reused salt.
    let (mut k_prime, theta) = hash::g(&ek_hash, &m_bytes, &salt);
    let (u2, v2) = pke::encrypt::<P>(&dk.ek, &m_bytes, &theta[..]);
    let c_prime = parsing::pack_ciphertext::<P>(&u2, &v2, &salt);

    // valid ⇔ decode succeeded AND the re-encryption reproduces c (CT compare).
    let reencrypt_ok = c_prime.as_slice().ct_eq(c);
    let valid = decode_ok & reencrypt_ok;

    // The selected key is the caller's to own; wipe the two stack key
    // candidates before returning. (Zeroize audit G2.)
    let shared = ct_select_key(&k_prime, &k_bar, valid);
    k_prime.zeroize();
    k_bar.zeroize();
    shared
}

// ── Helpers ───────────────────────────────────────────────────────────────────

/// Read `out.len()` bytes from a SHAKE256 reader. Thin wrapper so the XofReader
/// trait import stays local to this module.
fn read_xof(xof: &mut impl sha3::digest::XofReader, out: &mut [u8]) {
    xof.read(out);
}

/// Constant-time select: returns `a` when `choice` is set, else `b`. No branch
/// depends on `choice`.
fn ct_select_key(a: &SharedKey, b: &SharedKey, choice: Choice) -> SharedKey {
    let mut out = [0u8; SHARED_KEY_BYTES];
    for i in 0..SHARED_KEY_BYTES {
        out[i] = u8::conditional_select(&b[i], &a[i], choice);
    }
    out
}

// ── Constant-time audit: asm spot-check shim (Layer 3) ─────────────────────────
//
// `#[no_mangle] #[inline(never)]` so `ct_select_key` keeps a standalone symbol in
// `cargo asm` output. Compiled only under `--features ct-audit`; no effect on a
// normal build. The audit reads it to confirm the `subtle` select is branch-free
// masking (AND/XOR) and never a `jcc` on `choice`. See docs/audit/constant-time.md §5.
#[cfg(feature = "ct-audit")]
#[no_mangle]
#[inline(never)]
pub fn ct_asm_select_key(a: &SharedKey, b: &SharedKey, choice: u8) -> SharedKey {
    ct_select_key(a, b, Choice::from(choice & 1))
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::{Hqc128, Hqc192, Hqc256};

    // A tiny deterministic RNG so the random-entry-point paths are exercised
    // without pulling in an external rand implementation. NOT cryptographically
    // secure — test-only.
    struct TestRng(u64);
    impl TestRng {
        fn new(seed: u64) -> Self {
            TestRng(seed.wrapping_add(0x9E37_79B9_7F4A_7C15))
        }
    }
    impl RngCore for TestRng {
        fn next_u32(&mut self) -> u32 {
            (self.next_u64() >> 32) as u32
        }
        fn next_u64(&mut self) -> u64 {
            // SplitMix64.
            self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
            let mut z = self.0;
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            z ^ (z >> 31)
        }
        fn fill_bytes(&mut self, dest: &mut [u8]) {
            for chunk in dest.chunks_mut(8) {
                let r = self.next_u64().to_le_bytes();
                chunk.copy_from_slice(&r[..chunk.len()]);
            }
        }
        fn try_fill_bytes(&mut self, dest: &mut [u8]) -> Result<(), rand_core::Error> {
            self.fill_bytes(dest);
            Ok(())
        }
    }
    impl CryptoRng for TestRng {}

    // ── Core correctness: Decaps(Encaps(ek)) == K ─────────────────────────────

    fn kem_roundtrip<P: HqcParams>() {
        let seed_kem = [0x42u8; SEED_BYTES];
        let (ek, dk) = keygen_from_seed::<P>(&seed_kem);

        let m: Vec<u8> = (0..P::K).map(|i| (i as u8).wrapping_mul(7).wrapping_add(1)).collect();
        let salt = [0x17u8; SALT_BYTES];
        let (k_enc, c) = encaps_deterministic::<P>(&ek, &m, &salt);
        assert_eq!(c.len(), P::CT_BYTES);

        let k_dec = decaps::<P>(&dk, &c);
        assert_eq!(k_enc, k_dec, "shared keys must agree on the valid path");
    }

    #[test]
    fn kem_roundtrip_128() { kem_roundtrip::<Hqc128>(); }
    #[test]
    fn kem_roundtrip_192() { kem_roundtrip::<Hqc192>(); }
    #[test]
    fn kem_roundtrip_256() { kem_roundtrip::<Hqc256>(); }

    // ── Random entry points round-trip ────────────────────────────────────────

    #[test]
    fn kem_roundtrip_with_rng() {
        let mut rng = TestRng::new(1);
        let (ek, dk) = keygen::<Hqc128, _>(&mut rng);
        let (k_enc, c) = encaps::<Hqc128, _>(&mut rng, &ek);
        let k_dec = decaps::<Hqc128>(&dk, &c);
        assert_eq!(k_enc, k_dec);
    }

    // ── Determinism of keygen / encaps ────────────────────────────────────────

    #[test]
    fn keygen_from_seed_is_deterministic() {
        let seed = [0x5Au8; SEED_BYTES];
        let (ek1, dk1) = keygen_from_seed::<Hqc128>(&seed);
        let (ek2, dk2) = keygen_from_seed::<Hqc128>(&seed);
        assert_eq!(ek1.seed_ek, ek2.seed_ek);
        assert_eq!(ek1.s, ek2.s);
        assert_eq!(*dk1.sigma, *dk2.sigma);
        assert_eq!(*dk1.seed_kem, *dk2.seed_kem);
    }

    #[test]
    fn encaps_is_deterministic() {
        let (ek, _dk) = keygen_from_seed::<Hqc128>(&[1u8; SEED_BYTES]);
        let m = vec![3u8; Hqc128::K];
        let salt = [9u8; SALT_BYTES];
        let (k1, c1) = encaps_deterministic::<Hqc128>(&ek, &m, &salt);
        let (k2, c2) = encaps_deterministic::<Hqc128>(&ek, &m, &salt);
        assert_eq!(k1, k2);
        assert_eq!(c1, c2);
    }

    // ── Implicit rejection: a tampered ciphertext yields a different key ──────
    // It must not panic, must not return the original K, and must be stable
    // (the rejection key is deterministic in σ and c).

    fn implicit_rejection<P: HqcParams>() {
        let (ek, dk) = keygen_from_seed::<P>(&[0x33u8; SEED_BYTES]);
        let m: Vec<u8> = (0..P::K).map(|i| i as u8 ^ 0x5A).collect();
        let salt = [0x21u8; SALT_BYTES];
        let (k_valid, mut c) = encaps_deterministic::<P>(&ek, &m, &salt);

        // Flip one bit in the u-region (start of the ciphertext).
        c[0] ^= 0x01;

        let k_rej = decaps::<P>(&dk, &c);
        assert_ne!(k_valid, k_rej, "rejection key must differ from the valid K");

        // Deterministic: decapsulating the same tampered c again gives the same key.
        let k_rej2 = decaps::<P>(&dk, &c);
        assert_eq!(k_rej, k_rej2, "implicit-rejection key must be deterministic");

        // And it equals J(H(ek), σ, c) directly.
        let ek_hash = hash::h_ek(&ek.to_bytes());
        let expected = hash::j(&ek_hash, &dk.sigma, &c);
        assert_eq!(k_rej, expected, "rejection key must be J(H(ek), σ, c)");
    }

    #[test]
    fn implicit_rejection_128() { implicit_rejection::<Hqc128>(); }
    #[test]
    fn implicit_rejection_192() { implicit_rejection::<Hqc192>(); }
    #[test]
    fn implicit_rejection_256() { implicit_rejection::<Hqc256>(); }

    // ── Malformed ciphertexts never panic ─────────────────────────────────────

    #[test]
    fn decaps_handles_bad_length() {
        let (_ek, dk) = keygen_from_seed::<Hqc128>(&[7u8; SEED_BYTES]);
        // Empty, truncated, and oversized: each returns a key, no panic.
        let k0 = decaps::<Hqc128>(&dk, &[]);
        let k1 = decaps::<Hqc128>(&dk, &vec![0u8; Hqc128::CT_BYTES - 1]);
        let k2 = decaps::<Hqc128>(&dk, &vec![0u8; Hqc128::CT_BYTES + 1]);
        // They are J over the respective c, so empty vs truncated differ.
        assert_ne!(k0, k1);
        assert_ne!(k1, k2);
    }

    // ── Compressed secret key round-trip reproduces decaps behavior ───────────

    #[test]
    fn compressed_secret_key_roundtrip() {
        let seed = [0xC7u8; SEED_BYTES];
        let (ek, dk) = keygen_from_seed::<Hqc192>(&seed);
        let compressed = dk.to_bytes();
        assert_eq!(compressed.len(), SEED_BYTES);

        let dk2 = DecapsulationKey::<Hqc192>::from_bytes(&compressed[..]).expect("valid length");
        assert_eq!(*dk2.seed_kem, *dk.seed_kem);
        assert_eq!(*dk2.sigma, *dk.sigma);

        // A capsule made for the original key decapsulates identically under the
        // reconstructed key.
        let m = vec![0x11u8; Hqc192::K];
        let salt = [0x22u8; SALT_BYTES];
        let (k_enc, c) = encaps_deterministic::<Hqc192>(&ek, &m, &salt);
        assert_eq!(decaps::<Hqc192>(&dk, &c), k_enc);
        assert_eq!(decaps::<Hqc192>(&dk2, &c), k_enc);
    }

    #[test]
    fn from_bytes_rejects_bad_length() {
        assert!(DecapsulationKey::<Hqc128>::from_bytes(&[]).is_none());
        assert!(DecapsulationKey::<Hqc128>::from_bytes(&[0u8; SEED_BYTES - 1]).is_none());
    }
}
