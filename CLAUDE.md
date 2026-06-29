# CLAUDE.md — pqc-hqc

Context file for Claude Code. Read this before touching any code in this repo.

---

## Working instructions

- Never run `cargo` commands (e.g. `cargo test`, `cargo check`, `cargo build`). Only provide the command for the user to run in their terminal.
- Never run any code without asking. Usually it is better to provide/write some script and tell me to run it myself, giving me the command to do so.

---

## Project overview

Pure-Rust implementation of **HQC (Hamming Quasi-Cyclic)**, the NIST-selected
post-quantum KEM (selected March 2025, FIPS draft expected ~2027). Security reduces
to the Quasi-Cyclic Syndrome Decoding (QCSD) problem.

**Repo:** https://github.com/damesca/hqcr  
**License:** MIT OR Apache-2.0  
**Spec version:** [HQC specifications 22/08/2025](https://pqc-hqc.org/doc/hqc_specifications_2025_08_22.pdf)

---

## Directory layout

```
hqc/
├── src/
│   ├── lib.rs              # Crate root: re-exports, HqcParams trait
│   ├── params.rs           # Hqc128 / Hqc192 / Hqc256
│   ├── gf.rs               # GF(2^8) arithmetic (branch-free, table-free)
│   ├── poly/
│   │   ├── mod.rs          # Poly<P>: bit-packed [u64; N_WORDS]
│   │   ├── mul.rs          # Polynomial multiplication (hot path)
│   │   └── sampling.rs     # sample_fixed_weight / _mod / hamming_weight
│   ├── codes/
│   │   ├── mod.rs          # C.Encode / C.Decode (RS ∘ RM)
│   │   ├── reed_muller.rs  # RM(1,7): encode + duplicated FHT decode
│   │   └── reed_solomon.rs # Shortened RS: encode + CT BM + CT full-scan decode
│   ├── parsing.rs          # Serialize / deserialize keys and ciphertexts
│   ├── hash.rs             # G (SHAKE256), H (SHA3-256), J (SHAKE256), SHA3-512
│   ├── pke.rs              # HQC-PKE: Keygen, Encrypt, Decrypt (IND-CPA)
│   └── kem.rs              # HQC-KEM: Keygen, Encaps, Decaps (IND-CCA2, SFO⊥_m)
├── tests/
│   └── kat.rs              # KAT vectors — ground truth for correctness
├── benches/
│   └── bench.rs            # criterion: poly_mul, kem, decode_stages
└── docs/
    └── SECURITY-AUDIT.md   # CT + zeroize audit log
```

---

## Parameters (spec 2025 — authoritative)

> ⚠️ **`n1` = external RS code length**, **`n2` = internal RM code length**.
> Some earlier documents swap these — always check Table 5 of the spec.

| Constant | HQC-128 | HQC-192 | HQC-256 |
|:---|:---:|:---:|:---:|
| `n` (ring dimension) | 17 669 | 35 851 | 57 637 |
| `n1` (shortened RS length) | 46 | 56 | 90 |
| `n2` (duplicated RM length) | 384 | 640 | 640 |
| `k` (message bytes) | 16 | 24 | 32 |
| `ω` (secret weight) | 66 | 100 | 131 |
| `ωr = ωe` (ephemeral weight) | 75 | 114 | 149 |
| `δ` (RS error-correcting capacity) | 15 | 16 | 29 |
| RM multiplicity | 3 | 5 | 5 |
| `\|ekKEM\|` | 2 241 B | 4 514 B | 7 237 B |
| `\|cKEM\|` | 4 433 B | 8 978 B | 14 421 B |

---

## Mathematical core

### Ring and PKE

All arithmetic in `R = F2[X]/(X^n - 1)`, `n` primitive prime. Addition = XOR.

```
Secret key:  (x, y) ∈ R_ω × R_ω
Public key:  (h, s),  s = x + h·y

Encrypt:  u = r1 + h·r2,  v = C.Encode(m) + s·r2 + e
Decrypt:  v − u·y = C.Encode(m) + err;  m = C.Decode(·)
```

### RMRS codec

- **Encode:** RS(m → n1 symbols) then RM(each symbol → n2/n1 bits), concatenate → n1·n2 bits.
- **Decode:** RM decode each block (FHT) → n1 symbols; RS decode (BM + full scan) → m.
- RM base code: RM(1,7) = [128, 8, 64], duplicated ×3 or ×5.
- RS codes: `RS-S1: [46,16,31]`, `RS-S2: [56,24,33]`, `RS-S3: [90,32,49]`.
- Generator polynomials G1/G2/G3 are hardcoded from spec §3.4.2.

### KEM (SFO⊥_m)

```
G: SHAKE256 → (θ, K)    H: SHA3-256(ek)    J: SHAKE256 → K_bar

Decaps: valid = ct_eq(c'_PKE, c_PKE) & m'_not_none
        return ct_select(valid, K', J(σ, c))
```

---

## Key implementation decisions

### `poly/sampling.rs`

Two samplers:
- `sample_fixed_weight` — rejection sampling, for **x, y** (keygen only; non-CT, reference-accepted).
- `sample_fixed_weight_mod` — Barrett-reduction + backward swap, CT, for **r1, r2, e** (encrypt path).

### `poly/mul.rs`

- **Mode A** (`mul_sparse_dense`): sparse × dense, O(ω·N/64), used for keygen/encrypt.
- **Mode B** (`mul_dense_ct`): CT Karatsuba with `pclmulqdq` leaf (`clmul64`), used for decrypt (`u·y`). Portable fallback always present; SIMD enabled with `RUSTFLAGS="-C target-feature=+pclmulqdq"`.

### `gf.rs`

Branch-free, table-free: `gf_mul` uses 8-step carry-less multiply + `gf_reduce`; `gf_inv` uses addition chain `a^254`. `GF_EXP` retained for public-index paths only. `GF_LOG` unused (kept for potential future FFT).

### `codes/reed_solomon.rs`

CT decode pipeline: branchless syndromes → CT BM → `error_evaluator`/`formal_derivative` → full scan over all n1 positions (no Chien conditional push). Validity folded to one bit for the lone `Some`/`None` branch. All internal buffers `Zeroizing`.

---

## Constant-time requirements

| Location | What must be CT |
|:---|:---|
| `sampling.rs` `sample_fixed_weight_mod` | Fixed loops, `ct_eq` + `conditional_select` |
| `mul.rs` Mode B | Branchless bit extraction from secret `y` |
| `kem.rs::decaps` | `ct_eq(c'_PKE, c_PKE)`, `ct_select`, `Choice` for `m'_is_none` |
| `reed_muller.rs` argmax | 128-scan, `conditional_select` (no early exit) |
| `reed_solomon.rs` | Full-range scan, no syndrome fast path, no secret-indexed loads |

`sample_fixed_weight` (x/y) is **non-CT** — keygen-only, reference-accepted.

---

## Zeroize discipline

- `Poly<P>`: `#[derive(Zeroize, ZeroizeOnDrop)]`
- Seeds, σ, θ: `Zeroizing<[u8;N]>` wrappers
- `m'` in decaps: `Zeroizing<Vec<u8>>`
- RS/codec heap buffers: `Zeroizing` (syndromes, σ, Ω, σ', corrected, rs_cw, buf)
- RM decoder: stack arrays, not heap; no zeroize needed

---

## Dependencies

```toml
[dependencies]
sha3      = "0.10"
subtle    = "2.5"
zeroize   = { version = "1.7", features = ["derive"] }
rand_core = "0.6"

[dev-dependencies]
criterion = "0.5"
hex       = "0.4"
```

---

## Testing strategy

```bash
cargo test                        # unit + integration (fast)
cargo test --features kat         # KAT vectors — ground truth (~60s)
cargo test --features ct-audit --test ct_timing -- --nocapture
cargo bench
cargo bench --bench bench -- "kem/decaps" --baseline <name>
```

KAT tests (`tests/kat.rs`) are the **only** correctness oracle. If KAT passes for all three sets, the implementation is correct.

---

## Implementation status

All implementation steps complete. Current state:

| Step | Item | Status |
|:----:|:-----|:------:|
| 1–13 | Core modules (params → kem → lib) | ✅ done |
| 14 | Barrett `mod` sampler for r1/r2/e | ✅ done (KAT-verified) |
| 15 | KAT harness + vectors | ✅ done (byte-for-byte, all 3 sets) |
| 16 | `lib.rs` API polish | ✅ done |
| 17 | Karatsuba + `pclmulqdq` poly_mul | ✅ done |
| 18 | Criterion benchmarks | ✅ done |
| 19a/b | CT + zeroize audit | ✅ done (see `docs/SECURITY-AUDIT.md`) |
| 20a–e | CT RS/GF decoder (decoder `\|t\|` 622→0.82; G3 wiped) | ✅ done |
| 19c/d | Clippy/fmt, crate metadata | ⬜ deferred |
| 20c-opt | Additive-FFT root finder (perf-only) | ⬜ **not needed** — see perf note |
| 20-asm | Layer-3 asm/miri over Step-20 decoder | ⬜ follow-up |

### Performance (measured)

`pclmulqdq` restores decaps to near pre-CT levels. Without it the RS decoder dominates; with it `mul_dense_ct` dominates and the RS path is negligible. The decode_stages bench confirmed the FFT (20c-opt) would recover very little wall-clock time — **not worth implementing**.

| set | portable | +pclmulqdq | pre-CT master |
|---|---|---|---|
| hqc128 | ~2.46 ms | ~0.87 ms | ~0.78 ms |
| hqc192 | ~7.06 ms | ~2.22 ms | ~2.1 ms |
| hqc256 | ~17.07 ms | ~4.73 ms | ~4.4 ms |

Keygen and encaps are unaffected (they never decode).

---

## Reference implementations

- **Official C reference + AVX2:** https://pqc-hqc.org (source tarball)
- **pqcrypto-hqc:** Rust FFI bindings to C reference (not pure Rust)
- **RustCrypto/KEMs `hqc`:** placeholder `0.0.0`, not implemented

When spec text is ambiguous, the C reference is the tie-breaker.
