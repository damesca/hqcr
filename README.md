# hqc

A pure-Rust implementation of **HQC (Hamming Quasi-Cyclic)**, the code-based
post-quantum KEM selected by NIST in March 2025. Includes both the IND-CPA
public-key encryption scheme (HQC-PKE) and the IND-CCA2 key-encapsulation
mechanism (HQC-KEM) for all three parameter sets: **HQC-128 / HQC-192 / HQC-256**.

## ⚠️ Status — learning project, not production-ready

This crate was written as a **learning exercise** to understand the internals of
a modern code-based KEM end to end — finite-field arithmetic, concatenated
Reed-Muller / Reed-Solomon codes, quasi-cyclic polynomial multiplication, and
the salted Fujisaki–Okamoto transform.

**Do not use it to protect anything real.** In particular:

- It is **not validated against the official KAT (Known Answer Test) vectors**
  yet. The ephemeral sampler still uses rejection sampling instead of the spec's
  Barrett-reduction (`mod`) sampler, so the wire bytes do not yet match the
  reference implementation (see the roadmap in `CLAUDE.md`).
- It has **not** received an independent security review, side-channel audit, or
  formal constant-time verification.
- The implementation was developed **with the assistance of AI tools**, and has
  the limitations that implies — treat every line as needing review.

It *does* currently produce correct round-trips (encrypt→decrypt, encaps→decaps,
including implicit rejection) across all three parameter sets, verified by the
crate's unit tests.

## What is HQC?

HQC is a code-based IND-CCA2 KEM whose security reduces to the **Quasi-Cyclic
Syndrome Decoding** problem, a variant of the NP-complete Syndrome Decoding
problem. All arithmetic lives in the ring `R = F₂[X]/(Xⁿ − 1)` with `n` a
primitive prime. Messages are protected by a concatenated **Reed-Solomon ∘
Reed-Muller** code, and IND-CCA2 security comes from the salted FO transform
with implicit rejection (`SFO⊥_m`) layered over the IND-CPA PKE.

- **Spec:** [HQC specifications, 22/08/2025](https://pqc-hqc.org/doc/hqc_specifications_2025_08_22.pdf)
- **Paper:** Aguilar-Melchor et al., *Efficient Encryption from Random
  Quasi-Cyclic Codes*, IEEE Trans. Inf. Theory 64(5), 2018 (arXiv:1612.05572)

## Code structure

Single crate, bottom-up module layering. Each module depends only on the ones
above it.

```
src/
├── lib.rs              Crate root: module wiring + public re-exports
├── params.rs           HqcParams trait + Hqc128 / Hqc192 / Hqc256 constants
│
├── gf.rs               GF(2^8) arithmetic (log/exp tables) for Reed-Solomon
│
├── poly/
│   ├── mod.rs          Poly<P>: bit-packed ring element [u64; N_WORDS], XOR add
│   ├── mul.rs          Quasi-cyclic multiplication (sparse×dense, dense CT)
│   └── sampling.rs     Fixed-weight + uniform sampling from a SHAKE256 XOF
│
├── codes/
│   ├── mod.rs          C.Encode / C.Decode (the RS ∘ RM concatenation)
│   ├── reed_muller.rs  RM(1,7), duplicated; FHT-based decoder
│   └── reed_solomon.rs Shortened RS over GF(2^8); Berlekamp–Massey + Chien
│
├── parsing.rs          Wire-format (de)serialization of keys and ciphertexts
├── hash.rs             Keccak roles: G, H, J, the I seed-split, and the XOF
├── pke.rs              HQC-PKE: Keygen / Encrypt / Decrypt (IND-CPA)
└── kem.rs              HQC-KEM: Keygen / Encaps / Decaps (IND-CCA2, SFO⊥_m)
```

Design notes worth knowing:

- **Three parameter sets, always compiled.** Generic code is parameterized over
  the `HqcParams` trait; `Hqc128`, `Hqc192`, `Hqc256` are zero-sized type
  markers carrying the compile-time constants.
- **No heap for ring elements.** `Poly<P>` is a fixed `[u64; MAX_N_WORDS]` array
  on the stack, which keeps zeroization simple and avoids allocator noise.
- **Constant-time intent.** Secret-dependent paths (`kem::decaps` selection, the
  secret-times-ciphertext multiply, the code decoders) are written branchless
  using the `subtle` crate. This is *intent*, not *audited fact* — see the
  status note above.
- **Zeroization.** Secret key material, seeds, and ephemeral randomness are
  wrapped in `Zeroizing` / `ZeroizeOnDrop`.

For the full design rationale and the remaining roadmap (Barrett sampler → KAT
harness → API polish → SIMD `poly_mul`), see [`CLAUDE.md`](./CLAUDE.md).

## Usage

The KEM lives in the `hqc::kem` module and is generic over a parameter set.

### Encapsulate / decapsulate (with your own RNG)

The `keygen` / `encaps` entry points take any `rand_core::{RngCore + CryptoRng}`
source. The shared secret `K` produced by `encaps` matches the one recovered by
`decaps`.

```rust
use hqc::Hqc128;
use hqc::kem;

fn main() {
    // Bring your own cryptographically secure RNG (e.g. `rand::rngs::OsRng`,
    // which implements RngCore + CryptoRng).
    let mut rng = /* OsRng */ unimplemented!();

    // Keygen: public encapsulation key + secret decapsulation key.
    let (ek, dk) = kem::keygen::<Hqc128, _>(&mut rng);

    // Encapsulate against the public key → (shared secret, ciphertext).
    let (k_sender, ciphertext) = kem::encaps::<Hqc128, _>(&mut rng, &ek);

    // Decapsulate with the secret key → shared secret.
    let k_receiver = kem::decaps::<Hqc128>(&dk, &ciphertext);

    assert_eq!(k_sender, k_receiver); // both are [u8; 32]
}
```

### Deterministic API (reproducible, for tests / KAT)

Every operation also has a deterministic variant that takes the randomness
explicitly — handy for reproducibility and for the eventual KAT harness.

```rust
use hqc::{Hqc256, SEED_BYTES, SALT_BYTES};
use hqc::kem;

// Keygen from a fixed 32-byte seed_KEM.
let seed_kem = [0x42u8; SEED_BYTES];
let (ek, dk) = kem::keygen_from_seed::<Hqc256>(&seed_kem);

// Encaps from a fixed message (K bytes) and salt (16 bytes).
let m = [0x11u8; 32];            // Hqc256::K == 32
let salt = [0x22u8; SALT_BYTES];
let (k, ct) = kem::encaps_deterministic::<Hqc256>(&ek, &m, &salt);

assert_eq!(kem::decaps::<Hqc256>(&dk, &ct), k);
```

### Serialization

```rust
use hqc::Hqc192;
use hqc::kem;
use hqc::pke::EncryptionKey;

let seed = [7u8; hqc::SEED_BYTES];
let (ek, dk) = kem::keygen_from_seed::<Hqc192>(&seed);

// Public key ↔ bytes (ekKEM wire format: seed_ek ‖ s).
let ek_bytes = ek.to_bytes();
let ek_again = EncryptionKey::<Hqc192>::from_bytes(&ek_bytes).unwrap();

// Secret key is the compressed 32-byte seed_KEM; everything re-derives from it.
let dk_bytes = dk.to_bytes();                                  // Zeroizing<[u8; 32]>
let dk_again = kem::DecapsulationKey::<Hqc192>::from_bytes(&dk_bytes[..]).unwrap();
```

### Lower-level PKE

If you want the IND-CPA layer directly, `hqc::pke` exposes `keygen`, `encrypt`,
and `decrypt`. Note this layer is **not** CCA-secure on its own — use the KEM
for anything other than experimentation.

## Running the tests

The unit tests cover the GF arithmetic, both code layers, polynomial
multiplication, and the PKE / KEM round-trips (including implicit rejection and
malformed-ciphertext handling) for all three parameter sets.

```bash
cargo test
```

A KAT (Known Answer Test) feature is scaffolded but **not yet functional** — it
is blocked on the Barrett-reduction sampler. Once both land, the ground-truth
run will be:

```bash
cargo test --features kat   # not wired up yet — see CLAUDE.md roadmap
```

Benchmarks (criterion) are likewise scaffolded for a future step:

```bash
cargo bench                 # stub for now
```

## License

Licensed under either of **MIT** or **Apache-2.0** at your option.
