# hqcr

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

- It has **not** received an independent security review, side-channel audit, or
  formal constant-time verification.
- The implementation was developed **with the assistance of AI tools**, and has
  the limitations that implies — treat every line as needing review.

On **correctness**, though, it is in good shape. It is **validated against the
official NIST KAT (Known Answer Test) vectors**: the harness regenerates
`pk` / `sk` / `ct` / `ss` from the official `.req` seeds and asserts they match
the published `.rsp` files **byte-for-byte across all three parameter sets**. It
also reproduces the reference `intermediates_values` trace step-by-step, and the
encrypt→decrypt / encaps→decaps round-trips (including implicit rejection) pass
for every parameter set. See [Running the tests](#running-the-tests).

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
│   └── sampling.rs     Two fixed-weight samplers (rejection for x/y, Barrett
│                       "mod" for r1/r2/e) + uniform, from a SHAKE256 XOF
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

The dense×dense ring multiply (the decrypt hot path) uses limb-level
**Karatsuba** with a carry-less word-multiply leaf; on x86-64 built with
`+pclmulqdq` that leaf becomes a single `pclmulqdq` instruction (see
[Benchmarks](#benchmarks)). For the full design rationale and the remaining
roadmap (constant-time / zeroize audit, crate metadata), see
[`CLAUDE.md`](./CLAUDE.md).

## Building

### Default (portable, no `unsafe`)

The default build uses the portable Karatsuba + software carry-less multiply
leaf. Works on every target:

```bash
cargo build --release
cargo test
```

### With hardware `pclmulqdq` (x86-64, recommended on modern CPUs)

Enables the `pclmulqdq` carry-less multiply instruction as the Karatsuba leaf
— the only `unsafe` in the crate. Available on any x86-64 CPU since ~2010
(Intel Westmere / AMD Bulldozer and later).

```bash
# bash / Linux / macOS
RUSTFLAGS="-C target-feature=+pclmulqdq" cargo build --release

# PowerShell / Windows
$env:RUSTFLAGS="-C target-feature=+pclmulqdq"; cargo build --release
```

### Native CPU (all supported features, including `pclmulqdq`)

Lets the compiler use every instruction your current CPU supports. Produces a
binary that may not run on older machines:

```bash
# bash
RUSTFLAGS="-C target-cpu=native" cargo build --release

# PowerShell
$env:RUSTFLAGS="-C target-cpu=native"; cargo build --release
```

### KAT feature — byte-for-byte validation harness

The KAT and intermediate-values harnesses are gated behind the `kat` feature so
they don't ship in the library. Compile them in only for testing:

```bash
cargo test --features kat --test kat -- --nocapture
```

### Summary

| Build command | Karatsuba | Leaf multiply | `unsafe` |
|:---|:---:|:---:|:---:|
| `cargo build --release` | ✅ | portable (software) | none |
| `… -C target-feature=+pclmulqdq` | ✅ | `pclmulqdq` (hardware) | 1 block |
| `… -C target-cpu=native` | ✅ | `pclmulqdq` if available | 1 block |

## Usage

The KEM lives in the `hqcr::kem` module and is generic over a parameter set.

### Encapsulate / decapsulate (with your own RNG)

The `keygen` / `encaps` entry points take any `rand_core::{RngCore + CryptoRng}`
source. The shared secret `K` produced by `encaps` matches the one recovered by
`decaps`.

```rust
use hqcr::Hqc128;
use hqcr::kem;

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
explicitly — handy for reproducibility and used by the KAT harness.

```rust
use hqcr::{Hqc256, SEED_BYTES, SALT_BYTES};
use hqcr::kem;

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
use hqcr::Hqc192;
use hqcr::kem;
use hqcr::pke::EncryptionKey;

let seed = [7u8; hqcr::SEED_BYTES];
let (ek, dk) = kem::keygen_from_seed::<Hqc192>(&seed);

// Public key ↔ bytes (ekKEM wire format: seed_ek ‖ s).
let ek_bytes = ek.to_bytes();
let ek_again = EncryptionKey::<Hqc192>::from_bytes(&ek_bytes).unwrap();

// Secret key is the compressed 32-byte seed_KEM; everything re-derives from it.
let dk_bytes = dk.to_bytes();                                  // Zeroizing<[u8; 32]>
let dk_again = kem::DecapsulationKey::<Hqc192>::from_bytes(&dk_bytes[..]).unwrap();
```

### Lower-level PKE

If you want the IND-CPA layer directly, `hqcr::pke` exposes `keygen`, `encrypt`,
and `decrypt`. Note this layer is **not** CCA-secure on its own — use the KEM
for anything other than experimentation.

## Running the tests

### Unit tests

Cover the GF(2^8) arithmetic, both code layers (Reed-Muller, Reed-Solomon),
polynomial multiplication, the samplers, and the PKE / KEM round-trips
(including implicit rejection and malformed-ciphertext handling) for all three
parameter sets.

```bash
cargo test
```

The KAT and intermediate-values harnesses below are gated behind the `kat`
feature, so plain `cargo test` skips them.

### Known Answer Tests (KAT) — byte-for-byte validation

`tests/kat.rs` re-seeds the 2025 HQC KAT PRNG (a SHAKE256 XOF — *not* the legacy
AES-CTR DRBG) with each official `.req` seed, runs `keygen → encaps → decaps`,
writes a sibling `*.our.rsp`, and **asserts every `pk` / `sk` / `ct` / `ss`
matches the official `.rsp` byte-for-byte** for HQC-128 / 192 / 256. On a
mismatch it reports the count, field, and first differing offset.

```bash
cargo test --features kat --test kat -- --nocapture
```

The official vectors live in `tests/hqc-{1,3,5}/PQCkemKAT_*.rsp` and are the
assertion oracle — they are never overwritten.

### Intermediate-values trace

`tests/intermediate.rs` regenerates the reference `intermediates_values`
step-by-step trace (every keygen / encaps / decaps intermediate for count 0) and
writes `intermediates_values.our` alongside the official file.

```bash
cargo test --features kat --test intermediate -- --nocapture
```

It matches the reference `tests/hqc-{1,3,5}/intermediates_values` byte-for-byte;
to confirm, diff the two (they should be identical):

```bash
diff tests/hqc-1/intermediates_values tests/hqc-1/intermediates_values.our
```

### Benchmarks

Criterion benchmarks cover the ring-multiplication hot path (both modes) and the
full KEM operations, across all three parameter sets:

```bash
cargo bench                              # everything
cargo bench --bench bench -- poly_mul    # just the multiplication group
cargo bench --bench bench -- kem         # just keygen / encaps / decaps
```

The `poly_mul` group reports `sparse_dense` (Mode A, the keygen / encrypt path)
and `dense_ct` (Mode B, the constant-time decrypt path) per parameter set. The
`kem` group times `keygen` / `encaps` / `decaps`. HTML reports land in
`target/criterion/`.

Mode B uses limb-level Karatsuba with a carry-less word-multiply leaf. To
benchmark the hardware `pclmulqdq` variant (see [Building](#building)):

```bash
# bash
RUSTFLAGS="-C target-feature=+pclmulqdq" cargo bench --bench bench -- poly_mul

# PowerShell
$env:RUSTFLAGS="-C target-feature=+pclmulqdq"; cargo bench --bench bench -- poly_mul
```

Criterion compares against the previous run automatically, so run once without
the flag (portable baseline) then once with it to see the speedup on
`dense_ct/hqc128|192|256`.

## Sampler distribution analysis

An example binary generates empirical frequency distributions for the three
polynomial samplers and writes a self-contained HTML report you can open in any
browser.

```bash
cargo run --release --example sampler_distribution                 # HQC-128, 20 000 trials
cargo run --release --example sampler_distribution -- 192 50000    # HQC-192, 50 000 trials
cargo run --release --example sampler_distribution -- 256 20000    # HQC-256, 20 000 trials
```

Output: `sampler_distribution_<param>.html` in the current directory. Each
sampler card shows:

- **Per-position frequency curve** — should be flat at `weight/N` (green
  reference line); any positional bias appears as a slope or spike.
- **Histogram of per-position set-counts** — should match a
  `Binomial(trials, weight/N)` bell curve (orange overlay).
- **χ²/dof** — goodness-of-fit against the ideal flat distribution; values
  near 1.0 indicate an unbiased sampler.

`sample_fixed_weight` (secret x/y) and `sample_fixed_weight_mod` (ephemeral
r1/r2/e) should both read χ²/dof ≈ 1.0. The `mod` sampler may show a barely
visible excess at the very lowest positions (its backward duplicate-resolution
step maps collisions to small indices), but it is too small to see at typical
trial counts.

## License

Licensed under either of **MIT** or **Apache-2.0** at your option.
