# Fujisaki-Okamoto Transform

The Fujisaki-Okamoto (FO) transform upgrades an IND-CPA public-key encryption
scheme into an IND-CCA2 KEM. The core insight is that if the encryptor uses
*deterministic* randomness (derived from the message by a hash), the decapsulator
can **re-encrypt** and compare — turning any ciphertext tampering into a
detectable mismatch.

## Basic FO (1999 / 2013)

```
Encaps(ek):
    m  ←$ random
    θ   = H(m)              // derive randomness deterministically
    c   = PKE.Enc(ek, m; θ)
    K   = H'(m, c)
    return (K, c)

Decaps(dk, c):
    m'  = PKE.Dec(dk, c)
    θ'  = H(m')
    c'  = PKE.Enc(ek, m'; θ')
    if c' == c: return H'(m', c)
    else:       return ⊥
```

The re-encryption check is what makes the scheme CCA2-secure: an adversary
who flips bits in `c` causes `c' ≠ c`, so decapsulation fails cleanly instead
of leaking a partially-decrypted value.

**Cost:** one extra `PKE.Enc` inside every `Decaps` call.

## FO with Implicit Rejection (FO⊥)

Returning `⊥` on failure leaks a 1-bit oracle ("this ciphertext is invalid").
The fix — used in all NIST PQC KEMs — is *implicit rejection*: on failure,
return a pseudorandom key derived from a secret seed `σ` and the ciphertext
rather than a failure symbol.

```
Decaps(dk, c):
    m'  = PKE.Dec(dk, c)
    c'  = PKE.Enc(ek, m'; H(m'))
    if m' ≠ ⊥ AND c' == c: return KDF(m', c)
    else:                   return KDF(σ, c)   // implicit rejection
```

Both branches look identical to the caller: they always get a 32-byte key.
An adversary who submits a malformed ciphertext learns nothing, because
`KDF(σ, c)` is indistinguishable from a real key (σ is secret).

The comparison `c' == c` and the branch selecting `KDF(m', c)` vs `KDF(σ, c)`
must both be **constant-time** — the choice of which key is returned must not
be observable via timing.

## Salted FO with Implicit Rejection (SFO⊥_m) — 2025 HQC

The 2025 HQC spec uses a further hardened variant. Two additions:

**Salt.** A fresh 16-byte random salt is generated at encapsulation time and
included in the ciphertext. The salt is fed into `G` (the randomness-derivation
hash), so two encapsulations of the same message under the same key produce
different ciphertexts. This prevents multi-target and related-randomness attacks.

**Direct-message hashing.** The subscript `_m` means the message `m` is passed
directly to `G` (rather than being hashed first), keeping the PRF input uniform.

The resulting construction in HQC-KEM:

```
Encaps(ek):
    m    ←$ random(k bytes)
    salt ←$ random(16 bytes)
    ek_hash     = H(ek)                        // SHA3-256
    (θ, K)      = G(m, ek_hash, salt)          // SHAKE256
    (u, v)      = PKE.Enc(ek, m; θ)
    return (K, c = u || v || salt)

Decaps(dk, c):
    parse c = (u, v, salt)
    m'          = PKE.Dec(dk, (u, v))          // may be None
    ek_hash     = H(ek)
    (θ', K')    = G(m' or 0^k, ek_hash, salt) // always called — CT
    (u', v')    = PKE.Enc(ek, m'; θ')
    K_bar       = J(σ, c)                      // implicit-rejection key
    valid       = ct_not_none(m') & ct_eq(c_PKE', c_PKE)
    return ct_select(valid, K', K_bar)
```

`G`, `J` are SHAKE256 XOFs; `H` is SHA3-256. `ct_select` and `ct_eq` must
be constant-time (`subtle` crate in this implementation).

## Summary of variants

| Variant | Return on failure | Salt | Notes |
|:--------|:-----------------:|:----:|:------|
| FO (1999) | `⊥` (explicit) | no | original |
| FO⊥ | pseudorandom key | no | implicit rejection; used by Kyber/CRYSTALS |
| SFO⊥_m | pseudorandom key | yes | HQC 2025; salt in ciphertext |

The salt is the key difference between HQC's transform and Kyber's: it closes
a multi-ciphertext gap that is relevant for code-based schemes where the noise
distribution is less uniform than lattice schemes.
