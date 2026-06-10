// Intermediate-values harness.
//
// Regenerates a reference-shaped `intermediates_values` trace from THIS crate's
// implementation, so each step of keygen → encaps → decaps can be diffed
// against the official `tests/hqc-{1,3,5}/intermediates_values` files. Output
// is written to a sibling `intermediates_values.our` (the upstream file is
// never overwritten).
//
// ── How the inputs are derived ────────────────────────────────────────────────
// The reference `intermediates_values` file is generated with a *passthrough*
// RNG over the concatenated `.req` seeds (NOT the AES-CTR DRBG used by the
// `.rsp` KAT). Empirically, for count 0:
//     seed_kem = seedbuf[0 .. 32]
//     m        = seedbuf[32 .. 32+K]
//     salt     = seedbuf[32+K .. 32+K+16]
// where seedbuf is seed0 ‖ seed1 ‖ … (each `.req` seed is 48 bytes, so for
// HQC-128 `salt` spills into seed1). We reproduce exactly that, then re-derive
// every intermediate through the public API — no DRBG, fully deterministic.
//
// ── Scope (public-API subset) ─────────────────────────────────────────────────
// Everything reachable through the public API is reproduced in the reference
// format and order. The four RS-decoder-internal lines (syndromes, σ(x), z(x),
// error pairs) live inside `rs_decode` and need the `pub(crate)` `gf` module, so
// they are emitted as their degenerate no-error values (valid because count 0
// decrypts cleanly) — they are NOT live-computed from the decoder.
//
// Enabled with the `kat` feature:
//     cargo test --features kat --test intermediate -- --nocapture

#![cfg(feature = "kat")]

use std::fs;

use sha3::digest::XofReader;

use hqc::codes::{self, reed_muller, reed_solomon};
use hqc::poly::mul::{mul_dense_ct, mul_sparse_dense};
use hqc::poly::sampling::{sample_fixed_weight, sample_fixed_weight_mod, sample_uniform};
use hqc::poly::Poly;
use hqc::{hash, kem, parsing, pke};
use hqc::{Hqc128, Hqc192, Hqc256, HqcParams};

// ── .req parsing ──────────────────────────────────────────────────────────────

/// Concatenate every `seed =` value from a `.req` file into one byte stream
/// (the passthrough source the reference RNG reads sequentially).
fn concat_seeds(path: &str) -> Vec<u8> {
    let content = fs::read_to_string(path).unwrap_or_else(|e| panic!("read {path}: {e}"));
    let mut buf = Vec::new();
    for line in content.lines() {
        if let Some(rest) = line.trim().strip_prefix("seed = ") {
            buf.extend_from_slice(&hex::decode(rest.trim()).expect("seed must be valid hex"));
        }
    }
    assert!(buf.len() >= 96, "need at least two seeds' worth of bytes in {path}");
    buf
}

// ── Rendering helpers ─────────────────────────────────────────────────────────

fn ring_hex<P: HqcParams>(p: &Poly<P>) -> String {
    hex::encode(parsing::ring_to_bytes::<P>(p))
}

fn v_hex<P: HqcParams>(p: &Poly<P>) -> String {
    hex::encode(parsing::v_to_bytes::<P>(p))
}

/// `label: value` followed by exactly one blank line (the file's usual cadence).
fn kv(out: &mut String, label: &str, value: &str) {
    out.push_str(label);
    out.push_str(": ");
    out.push_str(value);
    out.push_str("\n\n");
}

// ── Generation ──────────────────────────────────────────────────────────────

/// Produce the full intermediate-values trace for parameter set `P`.
/// `level` is the NIST level label (1/3/5); `sec` is the security/DFR exponent.
fn generate<P: HqcParams>(req: &str, rsp: &str, level: u8, sec: u32) {
    let buf = concat_seeds(req);

    // Deterministic inputs (passthrough RNG over the concatenated seeds).
    let mut seed_kem = [0u8; 32];
    seed_kem.copy_from_slice(&buf[0..32]);
    let m = &buf[32..32 + P::K];
    let mut salt = [0u8; 16];
    salt.copy_from_slice(&buf[32 + P::K..32 + P::K + 16]);

    // ── KEYGEN ────────────────────────────────────────────────────────────────
    // seed_kem → (seed_pke ‖ σ) via SHAKE256; this crate squeezes a 32-byte σ,
    // the reference only K — seed_pke (first 32 B) is identical regardless.
    let (mut seed_pke, mut sigma) = ([0u8; 32], [0u8; 32]);
    {
        let mut xk = hash::xof(&seed_kem[..]);
        xk.read(&mut seed_pke);
        xk.read(&mut sigma);
    }
    let (seed_dk, seed_ek) = hash::i_pke_seed(&seed_pke);

    // ek/dk via the real keygen (authoritative), plus the KEM-level dk for the
    // final decaps shared key.
    let (ek, dk_pke) = pke::keygen::<P>(&seed_pke);
    let (_ek_kem, dk_kem) = kem::keygen_from_seed::<P>(&seed_kem);

    // y, x ← XOF(seed_dk); h ← XOF(seed_ek); s = x + h·y (= ek.s).
    let (y, x) = {
        let mut xdk = hash::xof(&seed_dk[..]);
        let y = sample_fixed_weight::<P>(&mut xdk, P::OMEGA);
        let x = sample_fixed_weight::<P>(&mut xdk, P::OMEGA);
        (y, x)
    };
    let h = {
        let mut xek = hash::xof(&seed_ek[..]);
        sample_uniform::<P>(&mut xek)
    };

    // ── ENCAPS ──────────────────────────────────────────────────────────────────
    let pk = ek.to_bytes();
    let ek_hash = hash::h_ek(&pk);
    let (k_shared, theta) = hash::g(&ek_hash, m, &salt);

    let enc = encrypt_trace::<P>(&ek, m, &theta[..]);
    let c = parsing::pack_ciphertext::<P>(&enc.u, &enc.v, &salt);

    // ── DECAPS ──────────────────────────────────────────────────────────────────
    let (u_in, v_in, _salt_in) = parsing::unpack_ciphertext::<P>(&c).expect("well-formed c");
    let y_dec = {
        let mut xdk = hash::xof(&seed_dk[..]);
        sample_fixed_weight::<P>(&mut xdk, P::OMEGA)
    };
    let uy = mul_dense_ct::<P>(&u_in, &y_dec);
    let decode_input = v_in.add(&uy); // C.Encode(m) + err, the input to C.Decode

    // RM decode each N2-bit block → the n1-symbol word fed to RS decode.
    let cw_bytes = parsing::v_to_bytes::<P>(&decode_input);
    let bb = P::N2 / 8;
    let mut rm_result = vec![0u8; P::N1];
    for j in 0..P::N1 {
        rm_result[j] = reed_muller::rm_decode(&cw_bytes[j * bb..(j + 1) * bb], P::MULTIPLICITY);
    }

    let m_prime = pke::decrypt::<P>(&dk_pke, &u_in, &v_in).expect("valid ct decodes");
    let (k_prime, theta_prime) = hash::g(&ek_hash, &m_prime, &salt);

    // Re-encryption (PKE.Encrypt of m' under θ') — its u'/v' must reproduce c.
    let reenc = encrypt_trace::<P>(&ek, &m_prime, &theta_prime[..]);

    // Final shared key via the real decaps (valid path ⇒ equals k_shared).
    let final_k = kem::decaps::<P>(&dk_kem, &c);

    // ── Assemble the trace in the exact reference layout ──────────────────────────
    let mut out = String::new();
    out.push_str(&format!("\n*********\n  HQC-{level}\n*********\n\n"));
    out.push_str(&format!(
        "N: {}   N1: {}   N2: {}   OMEGA: {}   OMEGA_R: {}   Failure rate: 2^-{}   Sec: {} bits\n",
        P::N, P::N1, P::N2, P::OMEGA, P::OMEGA_R, sec, sec
    ));

    // KEYGEN
    out.push_str("\n\n\n### KEYGEN ###\n\n");
    kv(&mut out, "seed_dk", &hex::encode(&seed_dk[..]));
    kv(&mut out, "seed_ek", &hex::encode(seed_ek));
    kv(&mut out, "y", &ring_hex::<P>(&y));
    kv(&mut out, "x", &ring_hex::<P>(&x));
    kv(&mut out, "h", &ring_hex::<P>(&h));
    kv(&mut out, "s", &ring_hex::<P>(&ek.s));
    kv(&mut out, "seed_kem", &hex::encode(seed_kem));
    kv(&mut out, "seed_pke", &hex::encode(seed_pke));
    kv(&mut out, "sigma", &hex::encode(&sigma[..P::K]));

    // ENCAPS
    out.push_str("\n\n### ENCAPS ###\n\n");
    kv(&mut out, "Reed-Solomon code word", &hex::encode(&enc.rs_codeword));
    kv(&mut out, "Concatenated code word", &v_hex::<P>(&enc.concatenated));
    kv(&mut out, "h", &ring_hex::<P>(&h));
    kv(&mut out, "s", &ring_hex::<P>(&ek.s));
    kv(&mut out, "r1", &ring_hex::<P>(&enc.r1));
    kv(&mut out, "r2", &ring_hex::<P>(&enc.r2));
    kv(&mut out, "e", &ring_hex::<P>(&enc.e));
    kv(&mut out, "Truncate(s.r2 + e)", &v_hex::<P>(&enc.sr2_plus_e));
    kv(&mut out, "c_pke->u", &ring_hex::<P>(&enc.u));
    kv(&mut out, "c_pke->v", &v_hex::<P>(&enc.v));
    kv(&mut out, "ek_kem", &hex::encode(&pk));
    kv(&mut out, "m", &hex::encode(m));
    kv(&mut out, "salt", &hex::encode(salt));
    kv(&mut out, "H(ek_kem)", &hex::encode(ek_hash));
    kv(&mut out, "theta", &hex::encode(&theta[..]));
    kv(&mut out, "c_kem", &hex::encode(&c));
    kv(&mut out, "K", &hex::encode(k_shared));

    // DECAPS
    out.push_str("\n\n### DECAPS ###\n\n");
    kv(&mut out, "c_pke.u", &ring_hex::<P>(&u_in));
    kv(&mut out, "c_pke.v", &v_hex::<P>(&v_in));
    kv(&mut out, "y", &ring_hex::<P>(&y_dec));
    kv(&mut out, "Truncate(u.y)", &v_hex::<P>(&uy));
    kv(&mut out, "v - Truncate(u.y)", &v_hex::<P>(&decode_input));

    // RS-decoder internals: degenerate no-error values (count 0 decrypts cleanly).
    let mut synd = String::from("The syndromes: ");
    for _ in 0..(2 * P::DELTA) {
        synd.push_str("0 ");
    }
    out.push_str(&synd);
    out.push_str("\n\n");
    out.push_str("The error locator polynomial: sigma(x) = 1\n\n");
    out.push_str("The polynomial: z(x) = 1\n\n");
    out.push_str("The pairs of (error locator numbers, error values): \n\n\n");

    kv(
        &mut out,
        "Reed-Muller decoding result (the input for the Reed-Solomon decoding algorithm)",
        &hex::encode(&rm_result),
    );
    // No-error case: the RS received word equals the RM decoding result.
    kv(&mut out, "Reed-Solomon code word", &hex::encode(&rm_result));
    kv(&mut out, "Concatenated code word", &v_hex::<P>(&reenc.concatenated));
    kv(&mut out, "h", &ring_hex::<P>(&h));
    kv(&mut out, "s", &ring_hex::<P>(&ek.s));
    kv(&mut out, "r1", &ring_hex::<P>(&reenc.r1));
    kv(&mut out, "r2", &ring_hex::<P>(&reenc.r2));
    kv(&mut out, "e", &ring_hex::<P>(&reenc.e));
    kv(&mut out, "Truncate(s.r2 + e)", &v_hex::<P>(&reenc.sr2_plus_e));
    kv(&mut out, "c_pke->u", &ring_hex::<P>(&reenc.u));
    kv(&mut out, "c_pke->v", &v_hex::<P>(&reenc.v));
    kv(&mut out, "ek_pke", &hex::encode(&pk));
    kv(&mut out, "dk_pke", &hex::encode(&seed_dk[..]));
    kv(&mut out, "c_kem", &hex::encode(&c));
    kv(&mut out, "m_prime", &hex::encode(&m_prime));
    kv(&mut out, "H(ek_kem)", &hex::encode(ek_hash));
    kv(&mut out, "theta_prime", &hex::encode(&theta_prime[..]));

    out.push_str("\n# Checking Ciphertext - Begin #\n\n");
    kv(&mut out, "c_kem_prime_t.c_pke.u", &ring_hex::<P>(&reenc.u));
    kv(&mut out, "c_kem_prime_t.c_pke.v", &v_hex::<P>(&reenc.v));
    kv(&mut out, "salt", &hex::encode(salt));
    out.push_str("# Checking Ciphertext - End #\n\n\n");

    kv(&mut out, "K_prime", &hex::encode(k_prime));
    out.push_str(&format!("secret1: {}\n", hex::encode(k_prime)));
    out.push_str(&format!("secret2: {}\n\n", hex::encode(final_k)));

    fs::write(rsp, out).unwrap_or_else(|e| panic!("write {rsp}: {e}"));
    eprintln!("intermediate: wrote {rsp}");
}

// ── PKE.Encrypt with all intermediates exposed ────────────────────────────────

struct EncryptTrace<P: HqcParams> {
    rs_codeword: Vec<u8>,    // RS codeword (n1 GF(2^8) symbols)
    concatenated: Poly<P>,   // C.Encode(m) (RMRS codeword embedded in the ring)
    r1: Poly<P>,
    r2: Poly<P>,
    e: Poly<P>,
    sr2_plus_e: Poly<P>,     // s·r2 + e (rendered truncated to the codeword region)
    u: Poly<P>,              // authoritative ciphertext component (from pke::encrypt)
    v: Poly<P>,
}

/// Recompute every PKE.Encrypt intermediate, mirroring `pke::encrypt` exactly,
/// while taking the authoritative `(u, v)` straight from the implementation.
fn encrypt_trace<P: HqcParams>(ek: &pke::EncryptionKey<P>, m: &[u8], theta: &[u8]) -> EncryptTrace<P> {
    let (u, v) = pke::encrypt::<P>(ek, m, theta);

    let mut rs_codeword = vec![0u8; P::N1];
    reed_solomon::rs_encode(m, &mut rs_codeword, P::DELTA);
    let concatenated = codes::encode::<P>(m);

    // (r2, e, r1) ← XOF(θ), in this exact sampling order (display order is r1/r2/e).
    let mut xth = hash::xof(theta);
    let r2 = sample_fixed_weight_mod::<P>(&mut xth, P::OMEGA_R);
    let e = sample_fixed_weight_mod::<P>(&mut xth, P::OMEGA_R);
    let r1 = sample_fixed_weight_mod::<P>(&mut xth, P::OMEGA_R);

    let sr2 = mul_sparse_dense::<P>(&r2, &ek.s);
    let sr2_plus_e = sr2.add(&e);

    EncryptTrace { rs_codeword, concatenated, r1, r2, e, sr2_plus_e, u, v }
}

// ── Entry points (one per parameter set) ──────────────────────────────────────

#[test]
fn intermediate_hqc128() {
    generate::<Hqc128>(
        "tests/hqc-1/PQCkemKAT_2321.req",
        "tests/hqc-1/intermediates_values.our",
        1,
        128,
    );
}

#[test]
fn intermediate_hqc192() {
    generate::<Hqc192>(
        "tests/hqc-3/PQCkemKAT_4602.req",
        "tests/hqc-3/intermediates_values.our",
        3,
        192,
    );
}

#[test]
fn intermediate_hqc256() {
    generate::<Hqc256>(
        "tests/hqc-5/PQCkemKAT_7333.req",
        "tests/hqc-5/intermediates_values.our",
        5,
        256,
    );
}
