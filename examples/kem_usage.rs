//! HQC-KEM usage example.
//!
//! Walks through the complete key-encapsulation flow for one or all three
//! parameter sets, covering:
//!
//!   1. Key generation, encapsulation, and decapsulation — the happy path.
//!   2. Implicit rejection: a tampered ciphertext yields a *different* but
//!      deterministic key; no panic, no error return.
//!   3. Compressed secret-key serialization / deserialization round-trip.
//!
//! Run (NEVER executed automatically — provided for the user to run):
//!
//!   cargo run --example kem_usage            # all three parameter sets
//!   cargo run --example kem_usage -- 128
//!   cargo run --example kem_usage -- 192
//!   cargo run --example kem_usage -- 256

use rand_core::OsRng;

use hqcr::kem::DecapsulationKey;
use hqcr::params::{Hqc128, Hqc192, Hqc256, HqcParams};
use hqcr::{decaps, encaps, keygen, SEED_BYTES, SHARED_KEY_BYTES};

fn main() {
    let arg = std::env::args().nth(1).unwrap_or_else(|| "all".to_string());
    match arg.as_str() {
        "128" => demo::<Hqc128>("HQC-128"),
        "192" => demo::<Hqc192>("HQC-192"),
        "256" => demo::<Hqc256>("HQC-256"),
        "all" => {
            demo::<Hqc128>("HQC-128");
            demo::<Hqc192>("HQC-192");
            demo::<Hqc256>("HQC-256");
        }
        other => {
            eprintln!("unknown parameter '{other}' — expected 128 | 192 | 256 | all");
            std::process::exit(1);
        }
    }
}

fn demo<P: HqcParams>(label: &str) {
    println!("=== {label} ===");
    println!("  Public key (ekKEM):      {:>6} bytes", P::PK_BYTES);
    println!("  Secret key (compressed): {:>6} bytes", SEED_BYTES);
    println!("  Ciphertext:              {:>6} bytes", P::CT_BYTES);
    println!("  Shared secret:           {:>6} bytes", SHARED_KEY_BYTES);
    println!();

    // ── Step 1: key generation ────────────────────────────────────────────────
    let (ek, dk) = keygen::<P, _>(&mut OsRng);

    // ── Step 2: encapsulation (sender) ────────────────────────────────────────
    // Returns a 32-byte shared secret K and the ciphertext to transmit.
    let (k_send, ct) = encaps::<P, _>(&mut OsRng, &ek);

    // ── Step 3: decapsulation (recipient) ─────────────────────────────────────
    // Always returns a 32-byte key. On success it equals k_send; on any
    // failure (bad ciphertext, decode error) it returns a deterministic
    // rejection key derived from the secret σ — indistinguishable by timing.
    let k_recv = decaps::<P>(&dk, &ct);

    assert_eq!(k_send, k_recv, "shared secrets must agree on the valid path");
    println!("  [1/3] keygen → encaps → decaps … OK");

    // ── Step 4: implicit rejection ────────────────────────────────────────────
    // Flip one byte anywhere in the ciphertext.
    let mut tampered = ct.clone();
    tampered[0] ^= 0xFF;
    let k_rej = decaps::<P>(&dk, &tampered);

    // The rejection key must differ from the valid key …
    assert_ne!(k_send, k_rej, "rejection key must differ from valid key");
    // … and must be deterministic: same tampered ciphertext → same key.
    assert_eq!(
        k_rej,
        decaps::<P>(&dk, &tampered),
        "implicit rejection must be deterministic"
    );
    println!("  [2/3] tampered ciphertext → implicit rejection … OK");

    // ── Step 5: compressed key round-trip ────────────────────────────────────
    // The full secret key is re-derived from a single 32-byte seed. Serialize
    // it, reconstruct the key, and confirm decapsulation still works.
    let seed_bytes = dk.to_bytes();
    assert_eq!(seed_bytes.len(), SEED_BYTES);
    let dk2 =
        DecapsulationKey::<P>::from_bytes(&seed_bytes[..]).expect("valid compressed secret key");
    assert_eq!(decaps::<P>(&dk2, &ct), k_send, "reconstructed key must agree");
    println!("  [3/3] compress → reconstruct → decaps … OK");

    println!();
}

