// KAT (Known Answer Test) harness.
//
// Reads the official NIST `.req` seed files from `tests/hqc-{1,3,5}/` and
// regenerates the `pk`, `sk`, `ct`, `ss` fields using THIS crate's
// implementation, writing the result to a sibling `*.our.rsp` file (the
// upstream `.rsp` vectors are never overwritten — they stay as the comparison
// oracle for manual inspection).
//
// ── How the seeds are consumed ────────────────────────────────────────────────
// Each `.req` `seed` is a 48-byte seed for the NIST AES-256-CTR DRBG
// (`randombytes`). We reproduce the standard NIST KAT call sequence per count:
//
//     randombytes_init(seed)
//     (pk, sk) = KEM.Keygen()      // draws SEED_BYTES from the DRBG
//     (ss, ct) = KEM.Encaps(pk)    // draws K + SALT_BYTES from the DRBG
//     ss'      = KEM.Decaps(sk, ct)
//     assert ss == ss'             // internal round-trip sanity check
//
// ── Caveat (read before judging the output) ───────────────────────────────────
// This crate's entropy layout differs from the HQC C reference, so the output
// will NOT match the upstream pqc-hqc.org `.rsp` byte-for-byte:
//   • Keygen here draws the 32-byte `seed_kem` from the DRBG, then expands it
//     via SHAKE256 into (seed_pke ‖ σ) and seed_pke into (seed_dk, seed_ek) —
//     the SAME architecture as the HQC reference (confirmed against the
//     `intermediates_values` trace). So given an identical `seed_kem` the keys
//     coincide; whether the DRBG hands us that same `seed_kem` is what the run
//     decides.
//   • σ length differs: this crate squeezes a 32-byte σ, the reference only K.
//     `seed_pke` (the first 32 bytes) is unaffected, so pk/seed_dk/seed_ek are
//     not; but the implicit-rejection path (which hashes σ) will diverge.
//
// The `sk` field is emitted in the reference wire format via
// `DecapsulationKey::expanded_secret_key_bytes()`:
//     sk = pk ‖ seed_dk (32) ‖ seed_ek (32) ‖ σ (K)        // 64 + K suffix
// matching the upstream |sk| (2321 / 4602 / 7333), NOT the compressed 32-byte
// `seed_kem`. The decaps round-trip (ss == ss') is asserted for every count, so
// the file is at minimum internally consistent.
//
// Enabled with the `kat` feature:
//     cargo test --features kat --test kat -- --nocapture

#![cfg(feature = "kat")]

use std::fs;

use aes::cipher::generic_array::GenericArray;
use aes::cipher::{BlockEncrypt, KeyInit};
use aes::Aes256;
use rand_core::{CryptoRng, RngCore};

use hqc::{kem, Hqc128, Hqc192, Hqc256, HqcParams};

// ── NIST AES-256-CTR DRBG (the `randombytes` used by the NIST KAT harness) ─────
//
// Faithful port of the reference `rng.c`: CTR_DRBG with AES-256 and no
// derivation function, 384-bit seed, single security strength. Each call to
// `randombytes` (here: each `fill_bytes`) emits the requested bytes from
// successive AES-CTR blocks and then runs the Update with no provided data.

struct NistDrbg {
    key: [u8; 32],
    v: [u8; 16],
}

impl NistDrbg {
    /// AES-256 single-block encryption (the DRBG's ECB primitive).
    fn aes256_ecb(key: &[u8; 32], input: &[u8; 16]) -> [u8; 16] {
        let cipher = Aes256::new(GenericArray::from_slice(key));
        let mut block = GenericArray::clone_from_slice(input);
        cipher.encrypt_block(&mut block);
        let mut out = [0u8; 16];
        out.copy_from_slice(block.as_slice());
        out
    }

    /// Big-endian increment of the 128-bit counter V (matches `rng.c`).
    fn increment_v(&mut self) {
        for j in (0..16).rev() {
            if self.v[j] == 0xff {
                self.v[j] = 0x00;
            } else {
                self.v[j] += 1;
                break;
            }
        }
    }

    /// AES256_CTR_DRBG_Update: refresh (Key, V), optionally XORing provided data.
    fn update(&mut self, provided: Option<&[u8; 48]>) {
        let mut temp = [0u8; 48];
        for i in 0..3 {
            self.increment_v();
            let block = Self::aes256_ecb(&self.key, &self.v);
            temp[16 * i..16 * i + 16].copy_from_slice(&block);
        }
        if let Some(pd) = provided {
            for i in 0..48 {
                temp[i] ^= pd[i];
            }
        }
        self.key.copy_from_slice(&temp[0..32]);
        self.v.copy_from_slice(&temp[32..48]);
    }

    /// randombytes_init(seed): zero (Key, V) then Update with the 48-byte seed.
    fn init(seed: &[u8; 48]) -> Self {
        let mut drbg = NistDrbg { key: [0u8; 32], v: [0u8; 16] };
        drbg.update(Some(seed));
        drbg
    }

    /// randombytes(out): fill `out` from AES-CTR blocks, then Update(None).
    fn randombytes(&mut self, out: &mut [u8]) {
        let mut i = 0;
        let mut remaining = out.len();
        while remaining > 0 {
            self.increment_v();
            let block = Self::aes256_ecb(&self.key, &self.v);
            if remaining > 15 {
                out[i..i + 16].copy_from_slice(&block);
                i += 16;
                remaining -= 16;
            } else {
                out[i..i + remaining].copy_from_slice(&block[..remaining]);
                remaining = 0;
            }
        }
        self.update(None);
    }
}

// Expose the DRBG as the RNG our KEM API consumes. Each `fill_bytes` is exactly
// one `randombytes` call (one Update at the end), matching the C reference's
// one-randombytes-per-buffer behaviour for keygen and encaps.
impl RngCore for NistDrbg {
    fn next_u32(&mut self) -> u32 {
        let mut b = [0u8; 4];
        self.randombytes(&mut b);
        u32::from_le_bytes(b)
    }
    fn next_u64(&mut self) -> u64 {
        let mut b = [0u8; 8];
        self.randombytes(&mut b);
        u64::from_le_bytes(b)
    }
    fn fill_bytes(&mut self, dest: &mut [u8]) {
        self.randombytes(dest);
    }
    fn try_fill_bytes(&mut self, dest: &mut [u8]) -> Result<(), rand_core::Error> {
        self.randombytes(dest);
        Ok(())
    }
}
impl CryptoRng for NistDrbg {}

// ── .req parsing ──────────────────────────────────────────────────────────────

/// Parse `(count, seed)` pairs from a NIST `.req` file. Only `count =` and
/// `seed =` lines carry data; `pk/sk/ct/ss` are blank in a `.req`.
fn parse_seeds(path: &str) -> Vec<(usize, [u8; 48])> {
    let content = fs::read_to_string(path).unwrap_or_else(|e| panic!("read {path}: {e}"));
    let mut out = Vec::new();
    let mut cur_count: Option<usize> = None;
    for line in content.lines() {
        let line = line.trim();
        if let Some(rest) = line.strip_prefix("count = ") {
            cur_count = Some(rest.trim().parse().expect("count must be an integer"));
        } else if let Some(rest) = line.strip_prefix("seed = ") {
            let bytes = hex::decode(rest.trim()).expect("seed must be valid hex");
            assert_eq!(bytes.len(), 48, "NIST DRBG seed must be 48 bytes");
            let mut s = [0u8; 48];
            s.copy_from_slice(&bytes);
            out.push((cur_count.expect("`count =` must precede `seed =`"), s));
        }
    }
    assert!(!out.is_empty(), "no seeds parsed from {path}");
    out
}

// ── Generation ──────────────────────────────────────────────────────────────

/// Run the full KAT sequence for parameter set `P` over every seed in `req` and
/// write the filled vectors to `rsp` (uppercase hex, NIST `.rsp` layout).
fn generate<P: HqcParams>(req: &str, rsp: &str, header: &str) {
    let seeds = parse_seeds(req);

    let mut output = String::new();
    output.push_str(&format!("# {header}\n\n"));

    for (count, seed) in &seeds {
        let mut drbg = NistDrbg::init(seed);

        // Standard NIST KAT call order: keypair, then enc, then dec.
        let (pk, dk) = kem::keygen::<P, _>(&mut drbg);
        let (ss, ct) = kem::encaps::<P, _>(&mut drbg, &pk);
        let ss_dec = kem::decaps::<P>(&dk, &ct);

        // Internal correctness oracle: the valid decapsulation must recover the
        // encapsulated key. This is what makes the generated file self-checking.
        assert_eq!(
            ss, ss_dec,
            "{header} count {count}: decaps(encaps()) shared-key mismatch"
        );

        let pk_bytes = pk.to_bytes();
        // Full reference-format sk (pk ‖ seed_dk ‖ seed_ek ‖ σ_K), not the
        // compressed 32-byte seed_KEM, so the field is comparable to upstream.
        let sk_bytes = dk.expanded_secret_key_bytes();
        assert_eq!(pk_bytes.len(), P::PK_BYTES, "pk length");
        assert_eq!(ct.len(), P::CT_BYTES, "ct length");

        output.push_str(&format!("count = {count}\n"));
        output.push_str(&format!("seed = {}\n", hex::encode_upper(seed)));
        output.push_str(&format!("pk = {}\n", hex::encode_upper(&pk_bytes)));
        output.push_str(&format!("sk = {}\n", hex::encode_upper(&sk_bytes[..])));
        output.push_str(&format!("ct = {}\n", hex::encode_upper(&ct)));
        output.push_str(&format!("ss = {}\n", hex::encode_upper(ss)));
        output.push('\n');
    }

    fs::write(rsp, output).unwrap_or_else(|e| panic!("write {rsp}: {e}"));
    eprintln!("KAT: wrote {rsp} ({} entries)", seeds.len());
}

// ── Entry points (one per parameter set) ──────────────────────────────────────

#[test]
fn kat_hqc128() {
    generate::<Hqc128>(
        "tests/hqc-1/PQCkemKAT_2321.req",
        "tests/hqc-1/PQCkemKAT_2321.our.rsp",
        "HQC-1",
    );
}

#[test]
fn kat_hqc192() {
    generate::<Hqc192>(
        "tests/hqc-3/PQCkemKAT_4602.req",
        "tests/hqc-3/PQCkemKAT_4602.our.rsp",
        "HQC-3",
    );
}

#[test]
fn kat_hqc256() {
    generate::<Hqc256>(
        "tests/hqc-5/PQCkemKAT_7333.req",
        "tests/hqc-5/PQCkemKAT_7333.our.rsp",
        "HQC-5",
    );
}
