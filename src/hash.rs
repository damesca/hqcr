// Hash and XOF wrappers around the sha3 crate (RustCrypto), instantiating the
// four Keccak roles of the HQC 2025 spec (SFO⊥_m transform) with their
// one-byte domain separators. Verified against the spec PDF (Table 1) and
// Saarinen's clean reference impl (github.com/mjosaarinen/hqc-py).
//
// Every Keccak call appends a single domain-separator byte AFTER the message
// (‖ = concatenation). The separators stop cross-function collisions in the
// random-oracle model — they are a security feature, not cosmetic.
//
//   XOF / seed expander :  SHAKE256( seed     ‖ 0x01 )      -> squeeze
//   H  (hash of ek)      :  SHA3-256( ek       ‖ 0x01 )      -> 32 B
//   I  (PKE seed split)  :  SHA3-512( seed_pke ‖ 0x02 )      -> 64 B
//                           split: seed_dk = [0:32], seed_ek = [32:64]
//   G  (K, θ derivation) :  SHA3-512( H(ek) ‖ m ‖ salt ‖ 0x00 ) -> 64 B
//                           split: K = [0:32], θ = [32:64]   (K FIRST, then θ)
//   J  (implicit reject) :  SHA3-256( H(ek) ‖ σ ‖ c ‖ 0x03 ) -> 32 B
//
// Note H and the XOF share the byte 0x01 but use different Keccak instances
// (SHA3-256 vs SHAKE256), so there is no collision between them.

use sha3::digest::{ExtendableOutput, Update, XofReader};
use sha3::{Digest, Sha3_256, Sha3_512, Shake256};
use zeroize::Zeroizing;

use crate::params::{SALT_BYTES, SEED_BYTES, SHARED_KEY_BYTES};

// ── Domain separators (spec Table 1) ──────────────────────────────────────────

const DOMAIN_XOF: u8 = 0x01; // SHAKE256 seed expander
const DOMAIN_H: u8 = 0x01; // SHA3-256 over ek
const DOMAIN_I: u8 = 0x02; // SHA3-512 PKE seed split
const DOMAIN_G: u8 = 0x00; // SHA3-512 (K, θ) derivation
const DOMAIN_J: u8 = 0x03; // SHA3-256 implicit-rejection key

// ── Sizes ───────────────────────────────────────────────────────────────────

/// Output length of H (SHA3-256) in bytes.
pub const DIGEST_BYTES: usize = 32;

/// Length of the ephemeral randomness seed θ produced by G, in bytes.
/// θ seeds the SHAKE256 XOF that samples (r2, e, r1) in PKE.Encrypt.
pub const THETA_BYTES: usize = 32;

/// 32-byte digest of the encapsulation key, H(ek). Public.
pub type EkHash = [u8; DIGEST_BYTES];

/// Ephemeral randomness seed θ. Secret — wrapped so it is zeroized on drop.
pub type Theta = Zeroizing<[u8; THETA_BYTES]>;

/// Final / candidate shared key K (also the implicit-rejection key). 32 bytes.
pub type SharedKey = [u8; SHARED_KEY_BYTES];

// ── Domain-separated Keccak helpers ────────────────────────────────────────────

/// SHA3-256 over `parts` concatenated, followed by the domain byte. 32-byte out.
fn sha3_256_ds(parts: &[&[u8]], domain: u8) -> [u8; 32] {
    let mut h = Sha3_256::default();
    for &p in parts {
        Digest::update(&mut h, p);
    }
    Digest::update(&mut h, [domain]);
    let digest = Digest::finalize(h);
    let mut out = [0u8; 32];
    out.copy_from_slice(&digest);
    out
}

/// SHA3-512 over `parts` concatenated, followed by the domain byte. 64-byte out.
/// The result is wrapped in `Zeroizing` because callers (G, I) split secret
/// material (θ, seed_dk) out of it.
fn sha3_512_ds(parts: &[&[u8]], domain: u8) -> Zeroizing<[u8; 64]> {
    let mut h = Sha3_512::default();
    for &p in parts {
        Digest::update(&mut h, p);
    }
    Digest::update(&mut h, [domain]);
    let digest = Digest::finalize(h);
    let mut out = Zeroizing::new([0u8; 64]);
    out.copy_from_slice(&digest);
    out
}

// ── H : SHA3-256(ek ‖ 0x01) ─────────────────────────────────────────────────────

/// H(ekKEM): the SHA3-256 digest of the encapsulation key. Computed once per
/// (de)capsulation and fed into G and J. Public input, no CT requirement.
pub fn h_ek(ek: &[u8]) -> EkHash {
    sha3_256_ds(&[ek], DOMAIN_H)
}

// ── G : SHA3-512(H(ek) ‖ m ‖ salt ‖ 0x00) → (K, θ) ──────────────────────────────

/// G derives the shared key K and the ephemeral randomness θ from the
/// encapsulation-key hash, the message, and the salt.
///
/// Returns `(K, θ)`: `K = digest[0:32]` (the shared key), `θ = digest[32:64]`
/// (the encryption randomness, kept secret and zeroized on drop).
pub fn g(ek_hash: &EkHash, m: &[u8], salt: &[u8; SALT_BYTES]) -> (SharedKey, Theta) {
    let digest = sha3_512_ds(&[ek_hash, m, salt], DOMAIN_G);

    let mut k = [0u8; SHARED_KEY_BYTES];
    k.copy_from_slice(&digest[..SHARED_KEY_BYTES]);

    let mut theta = Zeroizing::new([0u8; THETA_BYTES]);
    theta.copy_from_slice(&digest[SHARED_KEY_BYTES..SHARED_KEY_BYTES + THETA_BYTES]);

    (k, theta)
}

// ── J : SHA3-256(H(ek) ‖ σ ‖ c ‖ 0x03) → K_bar ──────────────────────────────────

/// J derives the implicit-rejection key from the encapsulation-key hash, the
/// secret σ, and the full ciphertext c. Returned by Decaps when the
/// re-encryption check fails; the caller selects between this and K' in
/// constant time.
pub fn j(ek_hash: &EkHash, sigma: &[u8; SEED_BYTES], c: &[u8]) -> SharedKey {
    sha3_256_ds(&[ek_hash, sigma, c], DOMAIN_J)
}

// ── XOF : SHAKE256(seed ‖ 0x01) ──────────────────────────────────────────────────

/// SHAKE256 XOF seeded by `seed ‖ 0x01`. Drives the fixed-weight and uniform
/// samplers: XOF(seed_ek) → h, XOF(seed_dk) → (y, x), XOF(θ) → (r2, e, r1), and
/// the KEM-level expansion XOF(seed_KEM) → (seed_pke, σ). The returned reader is
/// consumed by `poly::sampling::{sample_uniform, sample_fixed_weight}`.
pub fn xof(seed: &[u8]) -> impl XofReader {
    let mut x = Shake256::default();
    Update::update(&mut x, seed);
    Update::update(&mut x, &[DOMAIN_XOF]);
    x.finalize_xof()
}

// ── I : SHA3-512(seed_pke ‖ 0x02) → (seed_dk, seed_ek) ───────────────────────────

/// The PKE seed-split function I. Expands the 32-byte PKE seed into the dk-seed
/// (secret — feeds the sparse secret key (y, x), zeroized on drop) and the
/// ek-seed (public — feeds the uniform h). Called inside PKE.Keygen.
pub fn i_pke_seed(
    seed_pke: &[u8; SEED_BYTES],
) -> (Zeroizing<[u8; SEED_BYTES]>, [u8; SEED_BYTES]) {
    let digest = sha3_512_ds(&[seed_pke], DOMAIN_I);

    let mut seed_dk = Zeroizing::new([0u8; SEED_BYTES]);
    seed_dk.copy_from_slice(&digest[..SEED_BYTES]);

    let mut seed_ek = [0u8; SEED_BYTES];
    seed_ek.copy_from_slice(&digest[SEED_BYTES..2 * SEED_BYTES]);

    (seed_dk, seed_ek)
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use sha3::digest::{ExtendableOutput, Update, XofReader};
    use sha3::{Digest, Sha3_256, Sha3_512, Shake256};

    // ── Primitive wiring: NIST empty-input test vectors ───────────────────────
    // These pin that we are computing the correct Keccak instances, independent
    // of any HQC-specific concatenation or domain separation.

    #[test]
    fn sha3_256_empty_vector() {
        assert_eq!(
            hex::encode(Sha3_256::digest(b"")),
            "a7ffc6f8bf1ed76651c14756a061d662f580ff4de43b49fa82d80a4b80f8434a"
        );
    }

    #[test]
    fn sha3_512_empty_vector() {
        assert_eq!(
            hex::encode(Sha3_512::digest(b"")),
            "a69f73cca23a9ac5c8b567dc185a756e97c982164fe25859e0d1dcc1475c80a6\
             15b2123af1f5f94c11e3e9402c3ac558f500199d95b6d3e301758586281dcd26"
        );
    }

    // ── Domain separators are actually applied ────────────────────────────────

    #[test]
    fn xof_appends_domain_byte() {
        // xof(seed) must equal SHAKE256(seed ‖ 0x01), distinct from SHAKE256(seed).
        let seed = b"abc";
        let mut got = [0u8; 32];
        xof(seed).read(&mut got);

        let mut h = Shake256::default();
        Update::update(&mut h, seed);
        Update::update(&mut h, &[0x01u8]);
        let mut want = [0u8; 32];
        h.finalize_xof().read(&mut want);
        assert_eq!(got, want);

        // And it is NOT the plain (no-domain) SHAKE256.
        let mut h2 = Shake256::default();
        Update::update(&mut h2, seed);
        let mut plain = [0u8; 32];
        h2.finalize_xof().read(&mut plain);
        assert_ne!(got, plain);
    }

    #[test]
    fn h_ek_appends_domain_byte() {
        let ek = [0xABu8; 100];
        let mut buf = ek.to_vec();
        buf.push(0x01);
        assert_eq!(hex::encode(h_ek(&ek)), hex::encode(Sha3_256::digest(&buf)));
    }

    // ── G : SHA3-512(ek_hash ‖ m ‖ salt ‖ 0x00), split K=[0:32], θ=[32:64] ─────

    #[test]
    fn g_matches_spec_construction() {
        let ek_hash = [9u8; 32];
        let m = [7u8; 16];
        let salt = [5u8; 16];
        let (k, theta) = g(&ek_hash, &m, &salt);

        let mut buf = Vec::new();
        buf.extend_from_slice(&ek_hash);
        buf.extend_from_slice(&m);
        buf.extend_from_slice(&salt);
        buf.push(0x00);
        let expected = Sha3_512::digest(&buf);

        assert_eq!(&k[..], &expected[..32], "K must be the first 32 bytes");
        assert_eq!(theta.as_slice(), &expected[32..64], "θ must be the last 32 bytes");
    }

    #[test]
    fn g_deterministic_and_split_distinct() {
        let (k1, t1) = g(&[1u8; 32], &[2u8; 16], &[3u8; 16]);
        let (k2, t2) = g(&[1u8; 32], &[2u8; 16], &[3u8; 16]);
        assert_eq!(k1, k2);
        assert_eq!(*t1, *t2);
        assert_ne!(&k1[..], t1.as_slice());
    }

    #[test]
    fn g_sensitive_to_each_input() {
        let (bk, bt) = g(&[0u8; 32], &[0u8; 16], &[0u8; 16]);

        let mut ek = [0u8; 32];
        ek[0] = 1;
        let (k_ek, t_ek) = g(&ek, &[0u8; 16], &[0u8; 16]);

        let mut m = [0u8; 16];
        m[0] = 1;
        let (k_m, t_m) = g(&[0u8; 32], &m, &[0u8; 16]);

        let mut salt = [0u8; 16];
        salt[0] = 1;
        let (k_s, t_s) = g(&[0u8; 32], &[0u8; 16], &salt);

        for (k, t) in [(k_ek, t_ek), (k_m, t_m), (k_s, t_s)] {
            assert_ne!(bk, k);
            assert_ne!(*bt, *t);
        }
    }

    // ── J : SHA3-256(ek_hash ‖ σ ‖ c ‖ 0x03) ──────────────────────────────────

    #[test]
    fn j_matches_spec_construction() {
        let ek_hash = [0x44u8; 32];
        let sigma = [0x11u8; 32];
        let c = [0x22u8; 64];

        let mut buf = Vec::new();
        buf.extend_from_slice(&ek_hash);
        buf.extend_from_slice(&sigma);
        buf.extend_from_slice(&c);
        buf.push(0x03);
        assert_eq!(hex::encode(j(&ek_hash, &sigma, &c)), hex::encode(Sha3_256::digest(&buf)));
    }

    #[test]
    fn j_sensitive_to_each_input() {
        let ek_hash = [0x44u8; 32];
        let sigma = [0x11u8; 32];
        let c = [0x22u8; 64];
        let base = j(&ek_hash, &sigma, &c);

        let mut ek2 = ek_hash;
        ek2[0] ^= 1;
        assert_ne!(base, j(&ek2, &sigma, &c));

        let mut sigma2 = sigma;
        sigma2[0] ^= 1;
        assert_ne!(base, j(&ek_hash, &sigma2, &c));

        let mut c2 = c;
        c2[0] ^= 1;
        assert_ne!(base, j(&ek_hash, &sigma, &c2));
    }

    // ── I : SHA3-512(seed_pke ‖ 0x02), split seed_dk=[0:32], seed_ek=[32:64] ───

    #[test]
    fn i_pke_seed_matches_spec_construction() {
        let seed = [0x5Au8; SEED_BYTES];
        let (dk, ek) = i_pke_seed(&seed);

        let mut buf = seed.to_vec();
        buf.push(0x02);
        let full = Sha3_512::digest(&buf);
        assert_eq!(dk.as_slice(), &full[..SEED_BYTES]);
        assert_eq!(&ek[..], &full[SEED_BYTES..]);
        assert_ne!(dk.as_slice(), &ek[..]);
    }

    #[test]
    fn i_pke_seed_deterministic() {
        let seed = [0xC3u8; SEED_BYTES];
        let (dk1, ek1) = i_pke_seed(&seed);
        let (dk2, ek2) = i_pke_seed(&seed);
        assert_eq!(dk1.as_slice(), dk2.as_slice());
        assert_eq!(ek1, ek2);
    }
}
