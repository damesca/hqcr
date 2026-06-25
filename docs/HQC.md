# HQC — Hamming Quasi-Cyclic

> Succinct reference for this Rust implementation, based on the official
> specifications (v22/08/2025, `pqc-hqc.org`) and the foundational paper
> [AGM+2018]. Parameter names and sizes follow the spec convention used
> (**`n1` = outer RS length, `n2` = inner RM length**).

---

## 1. In one sentence

HQC is an ElGamal-like probabilistic encryption scheme over a ring of cyclic
binary polynomials, made IND-CCA2 via the Fujisaki-Okamoto transform. Its
security reduces to the **Quasi-Cyclic Syndrome Decoding (QCSD)** problem — a
structured variant of an NP-complete coding problem — and its decryption noise
is kept correctable by an error-correcting code.

Selected by NIST on **11 March 2025** as the code-based KEM to complement the
lattice-based ML-KEM (Kyber). Final FIPS expected ~2027.

---

## 2. The underlying algebra

### The ring R

All arithmetic lives in

$$\mathcal{R} = \mathbb{F}_2[X]/(X^n - 1)$$

where `n` is a **primitive prime** — chosen so that $(X^n-1)/(X-1)$ is
irreducible over $\mathbb{F}_2$. Hence $X^n - 1$ has exactly **two** irreducible
factors, which blocks algebraic attacks that would exploit a richer
factorization.

- Elements are degree-`< n` binary polynomials, i.e. vectors of $\mathbb{F}_2^n$
  (bit-packed as `[u64; N_WORDS]` here).
- **Addition** = bitwise XOR (no carries, no reduction).
- **Multiplication** = cyclic convolution mod $X^n - 1$:
  $w_k = \sum_{i+j \equiv k \,(\mathrm{mod}\ n)} u_i v_j$.

Multiplication by a fixed `v` equals a matrix-vector product with the
**circulant matrix** `rot(v)`. This is why HQC keys are compact: a quasi-cyclic
matrix is fully described by a single ring element of `n` bits, instead of the
$O(n^2)$ bits a random parity-check matrix would need.

### Fixed-weight vectors

$\mathcal{R}_\omega = \{v \in \mathcal{R} : \mathrm{wt}(v) = \omega\}$, the
sparse vectors of exact Hamming weight `ω`. The product of two sparse vectors
has bounded weight ($\mathrm{wt}(a\cdot b) \le \mathrm{wt}(a)\,\mathrm{wt}(b)$),
which is the lever that keeps decryption noise correctable.

---

## 3. The hard problem: QCSD

### Syndrome Decoding (the root)

Given a random $H \in \mathbb{F}_2^{(n-k)\times n}$ and a syndrome `y`, find a
low-weight `x` with $Hx^\top = y^\top$. This is **NP-complete**
(Berlekamp–McEliece–van Tilborg, 1978) and the best known attacks are
**Information Set Decoding (ISD)** variants (Prange → Stern/Dumer → MMT → BJMM),
exponential in `n`.

### The quasi-cyclic variant HQC uses

HQC relies on **2-DQCSD-P** (decisional, with a parity condition; HQC-256 adds
truncation → **3-DQCSD-PT**): distinguish $(H, y)$ where `H` is quasi-cyclic
(`rot`-structured) and `x` is sparse, from uniform. The QC structure shrinks
keys but only gives the attacker a sublinear advantage, which the parameter
sizing absorbs with margin. The parity condition kills trivial syndrome-parity
distinguishers at a cost of ≤ 1 bit.

---

## 4. HQC-PKE (IND-CPA)

ElGamal-like, with syndrome decoding playing the role of the trapdoor.

```
Secret key:  (x, y) ∈ R_ω × R_ω         // sparse, exact weight ω
Public key:  (h, s),  s = x + h·y,  h uniform in R

Encrypt(m; θ):                          // (r1, r2, e) ← R_ω sampled from θ
    u = r1 + h·r2                       // masked session value
    v = C.Encode(m) + s·r2 + e          // masked message

Decrypt(u, v):
    v - u·y = C.Encode(m) + (x·r2 + r1·y + e)
            = C.Encode(m) + err
    m = C.Decode(C.Encode(m) + err)     // succeeds iff wt(err) ≤ correction cap
```

By the public-key relation, $(h, s)$ is indistinguishable from uniform under
2-DQCSD-P. Decryption succeeds only when the noise
$err = x\cdot r_2 + r_1\cdot y + e$ stays within the code's correction
capacity; otherwise it fails. The **Decryption Failure Rate (DFR)** is required
$\le 2^{-\lambda}$ — both for correctness and because a high DFR is itself an
attack vector (chosen-ciphertext / reaction attacks on the secret key).

---

## 5. The RMRS concatenated code

The distinctive engineering piece of HQC: a **concatenated Reed-Muller ∘
Reed-Solomon** code recovers `m'` from `m' + err`.

- **Outer — shortened Reed-Solomon over GF(2⁸)** (`[46/56/90, k, …]`): corrects
  up to `δ` *symbol* errors via Berlekamp-Massey + Chien search. Uses the
  primitive polynomial `0x11D`.
- **Inner — Reed-Muller RM(1,7) = [128, 8, 64]**, duplicated 3× or 5× to
  [384] / [640]: each 8-bit RS symbol → one RM codeword; decoded with the
  **Fast Hadamard Transform** (argmax of the Walsh-Hadamard spectrum).

```
Encode:  m → RS.Encode → n1 symbols → RM.Encode each → concat (n1·n2 bits)
                                                        embedded in R, last ℓ bits truncated
Decode:  split into n1 blocks → FHT-decode each (RM) → RS.Decode (corrects symbol errors)
```

Concatenation gives much stronger effective correction than either code alone:
RM cleans up bit-errors inside a symbol, and RS recovers whole symbols even when
some RM blocks decode wrong.

---

## 6. HQC-KEM (IND-CCA2)

The IND-CPA PKE is lifted to an IND-CCA2 KEM with the **salted FO transform
with implicit rejection (SFO⊥_m)** introduced in the 2025 spec.

```
Keygen:   seed_KEM ←$ ; (seed_dk, seed_ek) = SHA3-512(seed_KEM)
          (ek, dk_PKE) = PKE.Keygen(seed_dk, seed_ek)
          dk = (ek, dk_PKE, σ, seed_KEM)         // compressed form: just seed_KEM (32 B)

Encaps(ek):
          m ←$ (k bytes) ; salt ←$ (16 bytes)
          (θ, K) = G(m, H(ek), salt)             // G = SHAKE256, H = SHA3-256
          c_PKE = PKE.Encrypt(ek, m; θ)
          return (K, c = c_PKE || salt)

Decaps(dk, c):
          m'   = PKE.Decrypt(dk_PKE, c_PKE)       // may fail
          (θ',K') = G(m' or 0, H(ek), salt)       // always called (constant-time)
          c'   = PKE.Encrypt(ek, m'; θ')          // re-encryption check
          valid = (m' ≠ ⊥) AND (c' == c_PKE)      // constant-time
          return valid ? K' : J(σ, c)             // implicit rejection
```

Two particularities worth highlighting:

- **Re-encryption** makes the ciphertext deterministic given `(m, ek)`; any
  tampering yields `c' ≠ c_PKE` and is detected — see [`FO.md`](FO.md).
- **Implicit rejection**: on failure, return a pseudorandom key `J(σ, c)`
  instead of an error, so the decapsulation oracle never reveals whether
  decryption succeeded.
- **Salt** (`salt` fed into `G`) randomizes encapsulation against multi-target
  attacks — the addition that distinguishes HQC's SFO⊥_m from Kyber's FO⊥.

The decryption-side selection (`valid ? K' : J`), the ciphertext comparison,
and the RM argmax / RS Chien search must all be **constant-time** — secret
material flows through them.

---

## 7. Parameters (spec 2025 — repo-authoritative)

| | HQC-128 | HQC-192 | HQC-256 |
|:--|:--:|:--:|:--:|
| NIST level | 1 | 3 | 5 |
| `n` (primitive prime) | 17 669 | 35 851 | 57 637 |
| `n1` (RS length, outer) | 46 | 56 | 90 |
| `n2` (RM length, inner) | 384 | 640 | 640 |
| `k` (message bytes) | 16 | 24 | 32 |
| `ω` (secret weight) | 66 | 100 | 131 |
| `ωr = ωe` (ephemeral weight) | 75 | 114 | 149 |
| `δ` (RS correction) | 15 | 16 | 29 |
| RM multiplicity | 3 | 5 | 5 |
| DFR target | < 2⁻¹²⁸ | < 2⁻¹⁹² | < 2⁻²⁵⁶ |
| public key `|ek|` | 2 241 B | 4 514 B | 7 237 B |
| ciphertext `|c|` | 4 433 B | 8 978 B | 14 421 B |
| secret key `|dk|` (compressed) | 32 B | 32 B | 32 B |
| shared key `|K|` | 32 B | 32 B | 32 B |

The secret key is just a 32-byte seed: `(x, y)` are regenerated on demand.
Compared to ML-KEM, HQC has larger ciphertexts but a tiny private key and an
**orthogonal** security base (coding theory vs. lattices) — the reason NIST
selected it as a hedge.

---

## 8. Main particularities at a glance

- **Compact keys from quasi-cyclicity**: a full circulant matrix collapses to
  one `n`-bit ring element.
- **Noise budget, not a trapdoor**: correctness rests on bounding the weight of
  sparse-vector products, decoded by RMRS — there is no algebraic trapdoor.
- **Two distinct samplers**: rejection sampling for the secret `(x, y)`; a
  constant-time Barrett (`mod`) sampler for the ephemeral `(r1, r2, e)` — the
  latter exists precisely to close the 2022 timing-attack class.
- **DFR is a security parameter**, not just a quality metric.
- **Conservative, transparent security**: best attack is ISD on the QC instance
  plus a sublinear structural factor, absorbed by parameter sizing.

---

## 9. References

- **[AGM+2018]** Aguilar-Melchor, Blazy, Deneuville, Gaborit, Zémor. *Efficient
  Encryption from Random Quasi-Cyclic Codes.* IEEE Trans. Inf. Theory 64(5),
  2018. arXiv:1612.05572. — the foundational paper.
- **[HQC-Spec-2025]** *Hamming Quasi-Cyclic (HQC), specifications v22/08/2025.*
  https://pqc-hqc.org — authoritative spec (SFO⊥_m, current parameters & DFR).
- **[HHK17]** Hofheinz, Hövelmanns, Kiltz. *A Modular Analysis of the
  Fujisaki-Okamoto Transformation.* TCC 2017.
- **[FS09]** Finiasz, Sendrier. *Security Bounds for the Design of Code-Based
  Cryptosystems.* ASIACRYPT 2009. — ISD complexity used for parameter sizing.

See also [`FO.md`](FO.md) for the FO transform in detail.
