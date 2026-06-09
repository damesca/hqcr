# CLAUDE.md — pqc-hqc

Context file for Claude Code. Read this before touching any code in this repo.

---

## Working instructions

- Never run `cargo` commands (e.g. `cargo test`, `cargo check`, `cargo build`). Only provide the command for the user to run in their terminal.

---

---

## Project overview

Pure-Rust implementation of **HQC (Hamming Quasi-Cyclic)**, the NIST-selected
post-quantum KEM (selected March 2025, FIPS draft expected ~2027). HQC is a
code-based IND-CCA2 KEM whose security reduces to the Quasi-Cyclic Syndrome
Decoding (QCSD) problem — a variant of the NP-complete Syndrome Decoding
problem.

This is a **pure-Rust, no-unsafe-by-default** implementation targeting
production quality: correct, constant-time, and verified against the official
KAT (Known Answer Test) vectors.

**Repo:** https://github.com/nshkrdotcom/pqc-hqc  
**License:** MIT OR Apache-2.0  
**Spec version:** [HQC specifications 22/08/2025](https://pqc-hqc.org/doc/hqc_specifications_2025_08_22.pdf)  
**Foundational paper:** Aguilar-Melchor et al., *Efficient Encryption from Random Quasi-Cyclic Codes*, IEEE Trans. Inf. Theory 64(5), 2018. arXiv:1612.05572

---

## Directory layout

Single-crate design. No workspace split — the internal complexity doesn't
justify the cross-crate compilation overhead for a standalone crypto library.

```
hqc/
├── Cargo.toml
├── CLAUDE.md
├── src/
│   ├── lib.rs              # Crate root: re-exports, HqcParams trait
│   ├── params.rs           # Hqc128 / Hqc192 / Hqc256 (all three always compiled)
│   │
│   ├── poly/
│   │   ├── mod.rs          # Poly type: bit-packed [u64; N_WORDS], add (XOR), reduce
│   │   ├── mul.rs          # Polynomial multiplication — the hot path
│   │   └── sampling.rs     # SampleFixedWeightVect$ + hamming_weight
│   │
│   ├── gf.rs               # GF(2^8) arithmetic for Reed-Solomon (log/exp tables)
│   │
│   ├── codes/
│   │   ├── mod.rs          # Top-level C.Encode / C.Decode (RS ∘ RM)
│   │   ├── reed_muller.rs  # RM(1,7): encode + duplicated FHT decode
│   │   └── reed_solomon.rs # Shortened RS over GF(2^8): encode + BM decode
│   │
│   ├── parsing.rs          # Serialize / deserialize keys and ciphertexts
│   ├── hash.rs             # G (SHAKE256), H (SHA3-256 of ek), J (SHAKE256), SHA3-512
│   ├── pke.rs              # HQC-PKE: Keygen, Encrypt, Decrypt (IND-CPA)
│   └── kem.rs              # HQC-KEM: Keygen, Encaps, Decaps (IND-CCA2, SFO⊥_m)
│
├── tests/
│   └── kat.rs              # KAT vectors — ground truth for correctness
│
└── benches/
    └── bench.rs            # criterion: poly_mul, keygen, encaps, decaps
```

---

## Parameters (spec 2025 — authoritative)

> ⚠️ The naming convention follows the spec: **`n1` = external RS code length**,
> **`n2` = internal RM code length**. Some earlier documents and other
> implementations swap these labels — always check against Table 5 of the spec.

| Constant | HQC-128 | HQC-192 | HQC-256 |
|:---|:---:|:---:|:---:|
| NIST security level | 1 | 3 | 5 |
| `n` (ring dimension, primitive prime) | 17 669 | 35 851 | 57 637 |
| `n1` (shortened RS code length) | 46 | 56 | 90 |
| `n2` (duplicated RM code length) | 384 | 640 | 640 |
| `k` (message length in bytes) | 16 | 24 | 32 |
| `ω` (secret vector weight) | 66 | 100 | 131 |
| `ωr = ωe` (ephemeral/error weight) | 75 | 114 | 149 |
| `δ` (RS error-correcting capacity) | 15 | 16 | 29 |
| RM multiplicity | 3 | 5 | 5 |
| DFR target | < 2⁻¹²⁸ | < 2⁻¹⁹² | < 2⁻²⁵⁶ |
| `\|ekKEM\|` (public key) | 2 241 B | 4 514 B | 7 237 B |
| `\|dkKEM\|` compressed | 32 B | 32 B | 32 B |
| `\|cKEM\|` (ciphertext) | 4 433 B | 8 978 B | 14 421 B |
| `\|K\|` (shared key) | 32 B | 32 B | 32 B |

**Derived sizes:**
- `|seed| = 32 B`, `|salt| = 16 B`
- `|ekKEM| = |seed| + ⌈n/8⌉`
- `|cKEM| = ⌈n/8⌉ + ⌈(n1·n2)/8⌉ + |salt|`

---

## Mathematical core

### The ring R

All arithmetic lives in `R = F2[X]/(X^n - 1)` where `n` is a **primitive prime**
(ensures `X^n - 1` has exactly two irreducible factors over F2, blocking
algebraic attacks). Elements are represented as bit-packed `[u64; N_WORDS]`.

```
w_k = sum_{i+j ≡ k (mod n)} u_i * v_j     // cyclic convolution
```

Addition in R = bitwise XOR. No carry, no reduction needed.

### HQC-PKE structure

```
Secret key:  (x, y) ∈ R_ω × R_ω          // sparse, weight exactly ω
Public key:  (h, s) where s = x + h·y     // h uniform random in R

Ciphertext:  (u, v) where
    u = r1 + h·r2                         // r1, r2, e ∈ R_ωr
    v = m' + s·r2 + e                     // m' = C.Encode(m)

Decryption:
    v - u·y = m' + (x·r2 + r1·y + e)
            = m' + err
    m = C.Decode(m' + err)  succeeds iff  weight(err) ≤ Δ
```

### RMRS codec — concatenated Reed-Muller / Reed-Solomon

**Encoding** (outer-then-inner, i.e. RS then RM):
1. RS.Encode: `m` (k bytes) → `n1` symbols over GF(2^8)
2. RM.Encode: each GF(2^8) symbol → one RM codeword of `n2/n1` bits (128 bits
   base, duplicated 3× or 5×)
3. Concatenate: total length = `n1 · (n2/n1)` = `n2` bits → embedded in R of
   length `n` (last `ℓ = n - n1·(n2/n1)` bits truncated)

> **Clarification on the codec parameters:** the RM base code is always
> **RM(1,7) = [128, 8, 64]** (m=7, blocks of 128 bits, encodes 8 bits per
> block). It is duplicated to [384, 8, 192] (×3) or [640, 8, 320] (×5).
> The RS code operates over **GF(2^8)** (m=8 for field arithmetic), using the
> primitive polynomial `1 + α² + α³ + α⁴ + α⁸`. These are two different uses
> of the letter `m` — do not conflate them.

**Decoding** (inner-then-outer, i.e. RM then RS):
1. Split into `n1` blocks of `(n2/n1)` bits each
2. RM.Decode each block via duplicated FHT → `n1` GF(2^8) symbols
3. RS.Decode via Berlekamp-Massey + Chien search → `m` (original message)

**Shortened RS codes** (Table 3 of spec):
```
RS-S1: [46,  16, 31]  (shortened from [255, 225, 31])
RS-S2: [56,  24, 33]  (shortened from [255, 223, 33])
RS-S3: [90,  32, 49]  (shortened from [255, 197, 49])
```
Precomputed generator polynomials `g1`, `g2`, `g3` are given in the spec
(§3.4.2) — hardcode them, do not recompute.

### HQC-KEM (IND-CCA2 via SFO⊥_m)

The transform used since the 2025 spec is the **salted FO with implicit
rejection**, denoted SFO⊥_m. Hash functions per the spec:

```
G : {0,1}* → {0,1}^|θ|×{0,1}^|K|     // SHAKE256, derives (θ, K) from (m, ek, salt)
H : {0,1}* → {0,1}^256                // SHA3-256, hashes ekKEM → used inside G
J : {0,1}* → {0,1}^|K|                // SHAKE256, rejection key from (σ, c)
```

Decaps always outputs a key — never an error:
```
if m' ≠ ⊥  AND  c'KEM == cKEM:   return K'        // valid path
else:                              return J(σ, c)   // implicit rejection
```

---

## Module responsibilities

### `params.rs`

Defines the `HqcParams` trait and three zero-size implementors:

```rust
pub trait HqcParams: sealed::Sealed {
    const N:           usize;  // ring dimension
    const N1:          usize;  // RS code length (external)
    const N2:          usize;  // RM code length after duplication (internal)
    const K:           usize;  // message length in bytes
    const OMEGA:       usize;  // secret weight
    const OMEGA_R:     usize;  // ephemeral weight ωr = ωe
    const DELTA:       usize;  // RS error-correcting capacity
    const MULTIPLICITY: usize; // RM duplication factor (3 or 5)
    // Derived:
    const N_WORDS:     usize;  // ceil(N / 64)
    const PK_BYTES:    usize;
    const CT_BYTES:    usize;
    const SK_BYTES:    usize;  // compressed = 32 B
}

pub struct Hqc128;
pub struct Hqc192;
pub struct Hqc256;
```

### `poly/mod.rs`

`Poly<P: HqcParams>`: bit-packed polynomial in R, stored as `[u64; P::N_WORDS]`
(stack-allocated, no heap). Operations:
- `add(&self, rhs) -> Self` — XOR (in-place variant too)
- `reduce(&mut self)` — fold the overflow bits back (only needed after sparse×dense mul)
- `get_bit(i)`, `set_bit(i)`, `clear()`

Do not use `Vec<u64>` — these arrays fit on the stack for all three parameter
sets and avoiding heap allocation matters for performance and zeroize.

### `poly/mul.rs`

**The hot path.** Two multiplication modes, both must be correct:

**Mode A — Sparse × Dense** (used for keygen and encrypt with secret/ephemeral
vectors): one operand has weight `ω` or `ωr`. Iterate over the set bit
positions and XOR-rotate the dense operand. Cost: `O(ω · N/64)`.

```
for pos in sparse.set_positions() {
    result ^= dense.rotate_left(pos);  // cyclic rotation by pos
}
```

**Mode B — Dense × Dense** (used in decrypt: `u · y` where `u` is an
arbitrary ciphertext component and `y` is a sparse secret): same as Mode A
but the sparse operand is the secret, so must be **constant-time** (the
positions are secret). Use Mode A with the sparse operand but with branchless
bit extraction.

Optimization layers (implement in order, gated by `#[cfg(target_feature)]`):
1. **L0 portable:** word-level, correct baseline
2. **L1 Karatsuba:** ~2× over L0 for Mode B dense×dense
3. **L2 SIMD:** `pclmulqdq` via `std::arch::x86_64` when `target_feature = "pclmul"`

### `poly/sampling.rs`

`sample_fixed_weight<P>(xof: &mut Shake256, n: usize, weight: usize) -> Poly<P>`:
rejection sampling that generates exactly `weight` distinct positions in `[0,n)`.
The loop count is variable but positions are chosen uniformly.

**CT requirement:** the number of rejection iterations must not branch on the
*value* of sampled positions (only on the public condition `pos < n`). Use
constant-time comparison for deduplication.

`hamming_weight(poly) -> usize`: popcount over `u64` limbs using
`u64::count_ones()`.

`sample_uniform<P>(xof: &mut Shake256) -> Poly<P>`: fills all `N` bits from
the XOF — used for the public `h`.

### `gf.rs`

GF(2^8) with primitive polynomial `p(x) = x⁸ + x⁴ + x³ + x² + 1`
(hex `0x11D`). Precompute two tables of 256 `u8` entries:
- `GF_LOG[256]`: discrete log base α (undefined for 0, set to sentinel)
- `GF_EXP[512]`: antilog, doubled to avoid modular reduction in multiply

```rust
pub fn gf_mul(a: u8, b: u8) -> u8 {
    if a == 0 || b == 0 { return 0; }
    GF_EXP[(GF_LOG[a as usize] as usize + GF_LOG[b as usize] as usize) % 255]
}
pub fn gf_inv(a: u8) -> u8 { GF_EXP[255 - GF_LOG[a as usize] as usize] }
```

All RS operations go through `gf.rs`. Zero is a special case everywhere.

### `codes/reed_muller.rs`

RM(1,7) base code: `[128, 8, 64]`. Encoding is a matrix-vector product over F2
(or equivalently: the 8-bit input directly indexes the 128-bit row of the
Hadamard matrix). Codeword is then duplicated `P::MULTIPLICITY` times.

**Decoding (duplicated FHT):**
1. Reshape the received block into `multiplicity` sub-blocks of 128 bits each.
2. For each position `i`, compute `F(i) = Σ_j (-1)^{v_{j,i}}` (sum over duplicates).
3. Apply the length-128 Walsh-Hadamard Transform (WHT) to the vector `F`.
4. The decoded 8-bit symbol is `argmax |F̂|`. The sign of `F̂` at that position
   determines the MSB (the all-ones correction vector).

Tie-breaking rule: if two positions have equal `|F̂|`, take the one with the
**smallest** value in the low 7 bits (per spec §3.4.3).

WHT works on `i16` intermediate values. Range check: with multiplicity 5 and
128 positions, values fit in `i16` without overflow.

**CT requirement:** the argmax loop must not early-exit. Always scan all 128
positions and update the running max with a branchless conditional move.

### `codes/reed_solomon.rs`

Shortened RS over GF(2^8). Uses the precomputed generator polynomials from
§3.4.2 of the spec — hardcode `G1`, `G2`, `G3` as `&[u8]` constants,
do not derive them at runtime.

**Encoding:** systematic; the message occupies the first `k` symbol positions
and parity fills the remaining `n1 - k` positions. Standard polynomial
long division over GF(2^8).

**Decoding:**
1. Compute `2δ` syndromes: `S_j = received(α^j)` for `j = 1..2δ`
2. Berlekamp-Massey → error locator polynomial `σ(x)`
3. Chien search: evaluate `σ(α^{-i})` for all `i ∈ [0, n1)` → error positions
4. Forney algorithm → error values (over GF(2^8); always 1 for binary codes
   but RS operates on 8-bit symbols so values are non-trivial)

**CT requirement:** Chien search must iterate over all `n1` candidates (46, 56,
or 90) without early exit, even after finding `deg(σ)` roots.

### `parsing.rs`

Byte-level serialization matching the spec wire format exactly (KAT vectors
validate this). Responsibilities:
- Bit-pack/unpack `Poly<P>` to/from `[u8; ceil(N/8)]`
- Pack public key: `seed_h (32 B) || s (ceil(n/8) B)`
- Pack ciphertext: `u (ceil(n/8) B) || v (ceil(n1·n2/8) B) || salt (16 B)`
- Compressed secret key: just `seed_KEM (32 B)`
- Full secret key: `ekKEM || dkPKE || σ (32 B) || seedKEM (32 B)`

Truncation: `v` covers exactly `n1 · n2` bits — where `n2` is the **duplicated**
RM length (384 for HQC-128, 640 for HQC-192/256), i.e. the full concatenated
RMRS codeword, *not* `n1 · 128`. This is embedded in the `n`-bit ring with the
trailing `ℓ = n - n1·n2` bits zeroed; those ℓ ring bits are dropped on
serialization. Cross-check: `CT_BYTES = ceil(n/8) + ceil(n1·n2/8) + 16` matches
the spec's published `|cKEM|` (4433 / 8978 / 14421) and `params.rs`.

### `hash.rs`

Thin wrappers around `sha3` (RustCrypto). All four roles from the spec:

```rust
// G: derives (θ, K) for encrypt and the final shared key
fn g(m: &[u8], ek_hash: &[u8; 32], salt: &[u8; 16]) -> (Theta, SharedKey)

// H: SHA3-256 hash of the encapsulation key (used inside G)  
fn h_ek(ek: &[u8]) -> [u8; 32]

// J: implicit rejection key
fn j(sigma: &[u8; 32], c: &[u8]) -> SharedKey

// Seed expander / XOF for sampling
fn xof(seed: &[u8]) -> impl XofReader   // SHAKE256
```

Do not use `SHA3-512` for key derivation — the spec 2025 uses `SHA3-256` for
`H(ekKEM)` and `SHAKE256` for `G` and `J`.

### `pke.rs`

```
HQC-PKE.Keygen(seed_dk, seed_ek):
    h ← sample_uniform(SHAKE256(seed_ek))
    (y, x) ← sample_fixed_weight × 2 (SHAKE256(seed_dk))
    s = x + h·y
    ek = (seed_ek, s),  dk = (seed_dk, ek)

HQC-PKE.Encrypt(ek, m, θ):
    (r2, e, r1) ← sample_fixed_weight × 3 (SHAKE256(θ))
    u = r1 + h·r2
    v = C.Encode(m) + s·r2 + e
    return (u, v)

HQC-PKE.Decrypt(dk, (u, v)):
    tmp = v + u·y             // tmp = C.Encode(m) + x·r2 + r1·y + e
    return C.Decode(tmp)      // returns None on decoding failure
```

Note on sampling order: the spec samples `(y, x)` (y first) from `seed_dk`
and `(r2, e, r1)` (r2 first) from `θ`. Match this order exactly — KAT vectors
will catch any swap.

### `kem.rs`

```
HQC-KEM.Keygen():
    seed_KEM ←$ random(32 B)
    σ ←$ random(32 B)
    (seed_dk, seed_ek) = SHA3-512(seed_KEM)
    (ek, dk_PKE) = PKE.Keygen(seed_dk, seed_ek)
    dk = (ek, dk_PKE, σ, seed_KEM)   // or compressed: (seed_KEM)

HQC-KEM.Encaps(ek):
    m ←$ random(k B)
    salt ←$ random(16 B)
    ek_hash = H(ek)
    (K, θ) = G(m, ek_hash, salt)
    c_PKE = PKE.Encrypt(ek, m, θ)
    return (K, c_PKE || salt)

HQC-KEM.Decaps(dk, c):
    parse c = (c_PKE, salt)
    m' = PKE.Decrypt(dk_PKE, c_PKE)    // may return None
    ek_hash = H(ek)
    (K', θ') = G(m' or zeros, ek_hash, salt)
    c'_PKE = PKE.Encrypt(ek, m', θ')
    K_bar = J(σ, c)
    // Constant-time select: valid iff m' ≠ None AND c'_PKE == c_PKE
    valid = (!m'_is_none) & ct_eq(c'_PKE, c_PKE)
    return ct_select(valid, K', K_bar)
```

When `m' = None` (decoding failure), use an all-zero buffer for the `G`
input — this ensures `G` is still called (constant-time) but `K_bar` will
be returned via `ct_select`.

---

## Constant-time requirements

Every location below is a security boundary. Non-CT code here is a
**vulnerability**, not a performance trade-off.

| Location | What must be CT | Use |
|:---|:---|:---|
| `poly/sampling.rs` | Rejection loop: position deduplication check | `subtle::ConstantTimeEq` for u16/u32 comparisons |
| `poly/mul.rs` Mode B | Bit extraction from secret sparse vector | Branchless `(word >> bit) & 1` over all positions |
| `kem.rs::decaps` | Ciphertext comparison `c'_PKE == c_PKE` | `subtle::ConstantTimeEq` byte-by-byte |
| `kem.rs::decaps` | Select `K'` vs `K_bar` | `subtle::ConditionallySelectable` |
| `kem.rs::decaps` | `m' == None` check feeds into select | `subtle::Choice`, not a bool branch |
| `codes/reed_muller.rs` | Argmax over 128 WHT outputs | Branchless running-max with `subtle::ConditionallySelectable` |
| `codes/reed_solomon.rs` | Chien search over all `n1` candidates | No early exit; always iterate full range |

Import `subtle::{ConstantTimeEq, ConditionallySelectable, Choice}`. Never use
`==` on anything derived from secret key material or decrypted message bytes.

---

## Zeroize discipline

Secret material that must be zeroized on drop:
- `Poly<P>` holding `x`, `y` (secret key components)
- `[u8; 32]` seeds
- `m` (recovered plaintext in Decaps)
- `θ` (ephemeral randomness)

Derive `ZeroizeOnDrop` from the `zeroize` crate on all types that hold the
above. Use `Zeroizing<[u8; N]>` wrappers for intermediate byte buffers.

---

## Dependencies

```toml
[dependencies]
sha3     = "0.10"        # SHAKE256, SHA3-256, SHA3-512
subtle   = "2.5"         # constant-time primitives
zeroize  = { version = "1.7", features = ["derive"] }
rand_core = "0.6"        # RngCore for seeding (optional, behind feature)

[dev-dependencies]
criterion = "0.5"        # benchmarks
hex       = "0.4"        # KAT vector parsing
```

No `unsafe` in any module except `poly/mul.rs` SIMD paths, gated behind
`#[cfg(target_feature = "...")]` with a safe portable fallback always present.

`std` target. A future `no_std + alloc` migration requires only renaming
`std::vec::Vec` → `alloc::vec::Vec` — no architectural changes needed.

---

## Testing strategy

```bash
cargo test                        # unit + integration tests (fast)
cargo test --features kat         # KAT vectors — slow (~60s), ground truth
cargo bench                       # criterion benchmarks
```

**KAT tests** (`tests/kat.rs`): parse official `.rsp` files from `pqc-hqc.org`
and verify `Keygen`, `Encaps`, `Decaps` output byte-for-byte for all three
parameter sets. These are the **only** correctness oracle that matters —
if KAT passes, the implementation is correct.

**Unit test coverage required:**
- `gf_mul` / `gf_inv` identities and commutativity
- `rm_encode` → inject up to `(192/2 - 1)` = 95 bit errors → `rm_decode` recovers
- `rs_encode` → inject up to `δ` symbol errors → `rs_decode` recovers
- `rs_decode` with `δ+1` errors → returns `None`
- `poly_mul` commutativity: `a·b == b·a` for 100 random pairs
- `pke_decrypt(pke_encrypt(m)) == m` for all three param sets
- `kem_decaps(kem_encaps(ek)) == K` for all three param sets
- Implicit rejection: flip one byte of a valid ciphertext → `decaps` returns
  a different key, does not panic, does not return the original `K`
- Decaps with zero-length or truncated ciphertext → no panic

---

## Implementation order

Bottom-up order based on module dependencies. Work through these one at a time.

| Step | File | Status | Depends on |
|:----:|:-----|:------:|:-----------|
| 1 | `src/params.rs` | ✅ done | — |
| 2 | `src/gf.rs` | ✅ done | — |
| 3 | `src/poly/mod.rs` | ✅ done | params |
| 4 | `src/poly/sampling.rs` | ✅ done | params, poly/mod |
| 5 | `src/poly/mul.rs` | ✅ done | params, poly/mod |
| 6 | `src/codes/reed_muller.rs` | ✅ done | params, poly/mod |
| 7 | `src/codes/reed_solomon.rs` | ✅ done | gf, params |
| 8 | `src/codes/mod.rs` | ✅ done | reed_muller, reed_solomon |
| 9 | `src/parsing.rs` | ✅ done | params, poly/mod |
| 10 | `src/hash.rs` | ✅ done | (sha3 crate only) |
| 11 | `src/pke.rs` | ✅ done | all of above |
| 12 | `src/kem.rs` | ✅ done | pke, hash |
| 13 | `src/lib.rs` | 🟡 module + re-exports wired | everything |

---

## Reference implementations

- **Official C reference + AVX2:** https://pqc-hqc.org (source tarball)  
  The AVX2 `poly_mul` shows the intended CLMUL/Karatsuba strategy.
- **pqcrypto-hqc:** Rust FFI bindings to the C reference (not pure Rust)
- **RustCrypto/KEMs `hqc`:** placeholder `0.0.0`, not implemented

When spec text is ambiguous, the C reference is the tie-breaker (after the
PDF spec). Sampling order and bit-packing conventions are easiest to verify
by running both implementations on the same seed and comparing intermediate
vectors.
