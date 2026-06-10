# KAT.md — Step 15 investigation notes

Status: **blocked, not done.** This file records everything discovered while
trying to make the official KAT vectors pass. It separates **verified facts**
(checked empirically against ground truth) from **open problems** (still
unresolved). The goal is to let a human verify the open items independently.

> Method note: I cannot run `cargo` in this environment, so I verified the
> algorithm logic in **Python** (`hashlib` has `sha3_256/sha3_512/shake_256`)
> against the official `tests/kat-vectors/hqc-128/intermediates_values` file,
> which prints the reference implementation's internal values for KAT case 0.
> Everything labelled "verified" below was reproduced byte-for-byte in Python.

---

## ★ Sampler fix (current session) — RESOLVED; SUPERSEDES §4 below

**Root cause (now confirmed against `reference_impl/src/vector.c`):** the
reference uses **two different fixed-weight samplers**, chosen per role in
`hqc.c`, and we had the secret-key one wrong.

| Reference routine | Used for | Algorithm |
|---|---|---|
| `vect_sample_fixed_weight1` (`…support1`) | **x, y** (keygen, hqc.c:47-48) | 24-bit **big-endian** draws in `3·weight`-byte batches; **reject** values ≥ `⌊2²⁴/N⌋·N`; reduce survivor `mod N` (Barrett); **redraw on duplicate** |
| `vect_sample_fixed_weight2` (`…support2`) | **r2, e, r1** (encrypt, hqc.c:113-115) | `4·weight` bytes as **little-endian u32**; `support[i] = i + ((rand·(N−i))>>32)`; backward dedup; no rejection |

Our `sample_fixed_weight_mod` was **already** an exact match for
`vect_sample_fixed_weight2`, so the ephemerals `r1/r2/e` were always going to be
fine. The bug was purely `x`/`y`: our old `sample_fixed_weight` was a bespoke
u16-LE rejection sampler that matched *neither* reference routine.

**Second root cause — `xof_get_bytes` is 8-byte aligned (symmetric.c):** every
XOF read rounds the request UP to a multiple of `sizeof(uint64_t)=8` bytes and
**discards** the unused tail of the final 8-byte unit:

```c
remainder = output_size % 8;
squeeze(output, output_size - remainder);          // aligned part
if (remainder) { squeeze(tmp, 8); copy `remainder` bytes; /* discard 8-remainder */ }
```

So each call advances the XOF by `ceil(output_size/8)·8` bytes. For `x`/`y` the
request is `3·66 = 198`; `198 % 8 = 6`, so each call consumes **200** bytes, not
198. `y` (first call) still matched because it only *uses* the first 198 — but
`x` then starts at offset **200**, not 198, so it was off by 2 bytes and wrong.
(Same for the ephemerals: `4·weight` is not a multiple of 8 for HQC-128/256, so
`e`/`r1` are misaligned after `r2`.) `h` was unaffected because its size
(`VEC_N_SIZE_64·8`) is already a multiple of 8 and it is the only draw on its XOF.

**Final change applied:**

| Site | Before (this session's wrong attempt) | After (correct) |
|---|---|---|
| `pke::keygen` — `y`, `x` | `sample_fixed_weight_mod` | `sample_fixed_weight` |
| `pke::decrypt` — re-derived `y` | `sample_fixed_weight_mod` | `sample_fixed_weight` |
| `src/poly/sampling.rs` — `sample_fixed_weight` body | u16-LE rejection | **rewritten** as reference `…support1` (24-bit BE, threshold reject, `% N`, redraw-on-dup, `3·weight` batched reads) |
| `src/poly/sampling.rs` — new `xof_get_bytes` helper | (reads were contiguous) | reads `len` bytes then discards `(8−len%8)%8` padding; used by **both** fixed-weight samplers |

(`r1/r2/e` stay on the mod sampler algorithm, now reading through `xof_get_bytes`
so `e`/`r1` align. `sample_uniform` for `h` is unchanged — its reads are already
8-aligned.)

**How the two-sampler split was pinned down (the investigation path):**
1. `tests/intermediate.rs` vs the reference `intermediates_values` showed
   `seed_dk`/`seed_ek`/`σ`/**`h`**/code-words all match, but `y`/`x` differ —
   isolating the bug to the fixed-weight sampler (the XOF stream is provably
   correct because `h`, read from the same XOF, matches).
2. `tests/sampler_probe.py` found the `y` position formula empirically:
   **24-bit big-endian `value % N`** (not 4-byte LE, not Barrett multiply-shift,
   no `+i`). This already proved *both* of our samplers were wrong for `x/y`.
3. `tests/sampler_probe2.py` confirmed the formula on `y` but showed `x` and the
   ephemerals misaligning — pointing at per-call rejection/redraw consuming a
   variable number of bytes.
4. Reading `reference_impl/src/vector.c` then made it exact: `…support1` is a
   **rejection** sampler (threshold + redraw-on-dup), and `…support2` is the mod
   sampler.
5. Implementing #1 fixed `y` but `x` was still wrong — which led to
   `reference_impl/src/symmetric.c`'s `xof_get_bytes`: it is 8-byte aligned and
   discards the tail, so calls after the first start at a padded offset. Adding
   that discard (the `xof_get_bytes` helper) is what aligns `x` and `e`/`r1`.

**Constant-time note:** `…support1` (hence our `sample_fixed_weight`) is
rejection-based and therefore **not constant time** — running time and XOF
consumption depend on the secret draws. This is inherent to matching the
reference (which accepts it for the once-per-keygen secret sampling); flagged for
the §19 CT audit.

> ⚠️ §4 below predates this and was left unresolved. It was *right* that neither
> of the then-existing samplers matched `y`; it just lacked the reference source
> to see there are two distinct routines. This section is authoritative.

---

## 0. The vector files (what each is)

For each parameter set (`hqc-128/192/256`) under `tests/kat-vectors/`:

- **`PQCkemKAT_*.rsp`** — the real KAT oracle. Per case: `seed` (48 B),
  `pk`, `sk`, `ct`, `ss`. The harness must reproduce these **byte-for-byte**.
- **`intermediates_values`** — HQC-specific debug dump for **case 0 only**:
  `seed_dk`, `seed_ek`, sampled vectors `x/y/h/s/r1/r2/e`, `m`, `salt`,
  `theta`, `K`, etc. Invaluable for locating *where* a mismatch is.
- `PQCkemKAT_*.req` — just the input seeds; redundant given the `.rsp`. (Not
  present after cleanup / can be ignored.)

⚠️ **Important subtlety discovered:** the `.rsp` and `intermediates_values`
were generated by **different harnesses** and do **not** share inputs (see §2).
The intermediates are good for verifying *math/derivation logic*; the `.rsp` is
the ground truth for the *full DRBG-driven flow*.

---

## 1. VERIFIED facts (reproduced in Python against intermediates, case 0)

All of the following match the reference **exactly**:

### 1.1 Hash / derivation chain
- `seed_kem → SHAKE256(seed_kem ‖ 0x01) → seed_pke (32 B) ‖ sigma (K B)`
  ✅ matches. **Note sigma length** — see open problem §3.1.
- `seed_pke → SHA3-512(seed_pke ‖ 0x02) → seed_dk (32 B) ‖ seed_ek (32 B)`
  ✅ matches our `hash::i_pke_seed` (domain `0x02`).
- `H(ek) = SHA3-256(ek ‖ 0x01)` ✅
- `G: SHA3-512(H(ek) ‖ m ‖ salt ‖ 0x00) → K (32 B) ‖ theta (32 B)`,
  split order **K first, theta second** ✅ matches our `hash::g`.

**Conclusion: `src/hash.rs` is correct.** No changes needed there except
nothing (the sigma length lives in `kem.rs`, not `hash.rs`).

### 1.2 Uniform sampler (public `h`)
- `h = sample_uniform(SHAKE256(seed_ek ‖ 0x01))` ✅ matches byte-for-byte.
- Confirms: the "seedexpander" for sampling **is** plain
  `SHAKE256(seed ‖ 0x01)`, read as a raw little-endian byte stream, with the
  bits above `N-1` masked off in the last word.
- **Conclusion: `sample_uniform` in `src/poly/sampling.rs` is correct, and the
  XOF domain byte `0x01` is correct for sampling too.**

### 1.3 NIST AES-256-CTR DRBG
- I implemented the NIST SP800-90A CTR_DRBG (no-df, AES-256), as used by the
  PQC `rng.c`, in pure Python and **self-tested it**: with the canonical fixed
  entropy `[0,1,2,…,47]`, `randombytes(48)` produces
  `061550234d158c5ec95595fe04ef7a25767f2e24cc2bc479d09d86dc9abcfde7056a8c266f9ef97ed08541dbd2e1ffa1`,
  which is the well-known NIST PQC first seed. ✅ DRBG implementation is correct.
- (The HQC `.rsp` seeds are *not* generated from `[0..47]` — HQC used a
  different fixed entropy to make the per-case seeds. That's fine; irrelevant
  to us. What matters is re-init'ing the DRBG with each case's 48-byte seed.)

---

## 2. OPEN PROBLEM — `.rsp` vs `intermediates` use different inputs

This caused a lot of confusion; documenting it so you don't fall into it.

For HQC-128 case 0, the `.rsp` `seed` is:
```
9EF877FDDBE8891C6E4E79EAF022E563DEFACA6B152161B9A423E8FE96A403E7  (bytes 0..31)
74B2D352CF74C934069C9DE74757F505                                  (bytes 32..47)
```

In **`intermediates_values`**:
- `seed_kem = 9ef877fd…96a403e7` = **`seed[0:32]` verbatim**
- `m        = 74b2d352…4757f505` = **`seed[32:48]` verbatim**
- (verified: `SHAKE256(seed_kem ‖ 0x01)` really does produce the
  intermediates' `seed_pke`/`sigma`, so `seed_kem = seed[0:32]` is genuinely
  the value the reference used here.)
- But `salt = aaf9baf4ae72c4c9b48efd574140a7bc` is **not** anywhere in the seed.

In **`.rsp`**:
- the case-0 shared secret `ss = 31D476B2A4D41B49…` **differs** from the
  intermediates' `K = a753321c2cb26174…`.

**Interpretation:** the `intermediates_values` harness sets
`seed_kem = seed[0:32]` and `m = seed[32:48]` *directly* (a debug shortcut),
whereas the `.rsp` harness uses the **standard NIST flow**: re-init the DRBG
with the 48-byte seed, then draw `seed_kem`, `m`, `salt` via `randombytes`.
I confirmed the DRBG output for the `.rsp` seed:
```
randombytes(32) after init(seed) = a01abd3292c0fbc23a39125421852a41204e891560c34178cbd0ebb7f575f31d
```
which is **not** equal to `seed_kem` in the intermediates (so the intermediates
did not use the DRBG). 

**What still needs confirming (the open part):** the exact `randombytes` call
**sequence and sizes** the `.rsp` harness used. The standard NIST HQC
`PQCgenKAT_kem.c` flow is believed to be:
```
randombytes_init(seed_48, NULL, 256)
crypto_kem_keypair:  randombytes(seed_kem, 32)         # 1 draw
crypto_kem_enc:      randombytes(m, K)                  # then
                     randombytes(salt, SALT_BYTES=16)   # (order m-then-salt assumed)
```
I could **not** verify this end-to-end because it requires running the full
keygen+encaps (sampling + poly mul + RS/RM encode), which I can't do in Python
quickly and can't run in cargo. **This ordering must be confirmed** (it changes
the DRBG state and therefore every output byte). The intermediates file does
not help here because it bypassed the DRBG.

---

## 3. OPEN PROBLEM — wire-format / size mismatches with the reference

### 3.1 `sigma` is **K bytes, not 32**
- Verified: the XOF split is `seed_pke (32) ‖ sigma (K)` where K = 16/24/32 for
  HQC-128/192/256. In case 0, `sigma = d2f55cdd7665ebf10318f96adc8dc30c`
  (16 bytes).
- Our `kem.rs` currently reads `sigma` as `[u8; SEED_BYTES]` (32 B) and reads it
  *after* `seed_pke` from the same XOF. For HQC-128 this over-reads by 16 bytes.
- Impact: only affects the implicit-rejection key `J`. The standard `.rsp`
  tests the **valid** decaps path (where `ss = K`, independent of sigma), so
  this may not break the `.rsp` `ss`/`ct`/`pk` — **but** it is wrong vs the
  reference and must be fixed for correctness. Fix: make `sigma` `K` bytes.

### 3.2 `.rsp` `sk` is the **full 2321-byte reference format**, not our 32 B
- Measured lengths (HQC-128 case 0): `pk = 2241 B`, `sk = 2321 B`,
  `ct = 4433 B`, `ss = 32 B`.
- `2321 = 2241 (pk) + 80`. The `.rsp` `sk` **begins with the full pk bytes**
  (`sk` and `pk` share the prefix `4053237912EA281C…`). So the reference `sk`
  layout is roughly `sk = <secret seed(s)> ‖ … ‖ pk` or `pk ‖ <secrets>` —
  **exact layout unconfirmed** (the 80 trailing/leading bytes need to be
  decoded: plausibly `sk_seed (32) ‖ sigma ‖ pk`, but `32 + 16 = 48 ≠ 80`, so
  the seed sizes used by the reference may differ — possibly a 40-byte
  seedexpander seed; **needs confirmation against the reference `kem.c`**).
- Impact: **our `DecapsulationKey::to_bytes()` returns only the compressed
  32-byte `seed_kem`.** To match the `.rsp` `sk` byte-for-byte we must add a
  KAT-format `sk` serializer that reproduces the reference layout exactly.
  (Alternatively, the KAT harness could skip the `sk` comparison and only check
  `pk`/`ct`/`ss` — but that is weaker.)

### 3.3 `pk` format
- `pk = 2241 B = seed_ek (32) ‖ s (2209)` is consistent with our
  `parsing::pack_public_key`. Not yet verified byte-for-byte against `.rsp`
  `pk` because that needs a full keygen run, but the structure matches and
  `seed_ek`/`s` were verified individually in the intermediates.

---

## 4. OPEN PROBLEM (CRITICAL) — fixed-weight sampler for `x`, `y` is wrong

This is the main blocker. **Neither** of our two samplers reproduces the
intermediates' `y` (or `x`):

- `sample_fixed_weight` (rejection, u16 LE, reject ≥ N) — ❌ no match
- `sample_fixed_weight_mod` (Barrett, u32 LE, `i + ((rand·(N−i))>>32)`,
  backward dedup) — ❌ no match

### What I tried (all against `y` and `x`, set-equality, case 0)
Stream confirmed correct: `SHAKE256(seed_dk ‖ 0x01)` (same construction that
made `h` match). Variants tested:
- read width **3 bytes** and **4 bytes** per candidate
- endianness **little** and **big**
- reduction **Barrett** (`i + ((rand·(N−i)) >> 8·width)`) and **true modulo**
  (`i + (rand mod (N−i))`)
- also `rand·N >> …` and `rand mod N` (no `+i` offset)
- dedup: none, and backward-collision→`i`
- stream offsets: `0` (y sampled first), `4·ω`, `3·ω` (in case x is first)

**None matched.** (The last comprehensive sweep was interrupted before
printing, but earlier targeted runs of the most likely candidates — Barrett
u32-LE, mod 3-byte-BE — definitively did not match `y`'s bit set.)

### Ground-truth data for `y` (HQC-128 case 0)
- `seed_dk = 12daf031bdc7fc592e0003a21eefa9a1019539abccc8f67075947cbfeaac98c5`
- `y` set bit positions (weight 66, sorted):
```
50, 123, 229, 414, 506, 765, 936, 1537, 1539, 1691, 1833, 1931, 1939, 1961,
2082, 2495, 2627, 2842, 3887, 4298, 4473, 4712, 4796, 4897, 5477, 5933, 5996,
6314, 6604, 6724, 7248, 7551, 7576, 8478, 8535, 8552, 8642, 8947, 9486, 9918,
9940, 10227, 10479, 11277, 11735, 12260, 12336, 12492, 12650, 13041, 13141,
13226, 13442, 13578, 13892, 13949, 14340, 14542, 14979, 15513, 15708, 15823,
16178, 16614, 16773, 17109
```
(Positions are sorted only because they are read out of a bit-vector; this is
**not** the order the sampler produced them in.)

### Hypotheses for why it doesn't match (to investigate)
1. **The reference fixed-weight sampler uses a different reduction or read
   width** than what we implemented. The exact `vect_set_random_fixed_weight`
   from the **2025 reference C source** (pqc-hqc.org tarball) must be read
   line-by-line. There are multiple historical variants; ours may be from the
   wrong revision.
2. **The seedexpander may buffer/emit in fixed blocks** so that the bytes
   consumed by the fixed-weight sampler are not a contiguous `stream[0:4ω]`
   slice. (Unlikely — `h` matched a contiguous slice — but the fixed-weight
   path requests a different number of bytes and the reference's seedexpander
   has a 64-byte internal buffer in some versions; worth checking.)
3. **Sampling order**: maybe `x` is sampled before `y`, so `y` uses a later
   slice of the stream. (Tested offsets `4ω`/`3ω`; still no match, but only for
   the variants tried.)
4. **Possible extra domain/diversifier** in the seedexpander init for the
   secret-key expander vs the public one. (`h` from `seed_ek` used `0x01`;
   maybe `seed_dk` path differs — though that would be unusual.)

### Recommended way to resolve (for the human)
Get the **2025 reference C** (`pqc-hqc.org` tarball, file
`src/vector.c` → `vect_set_random_fixed_weight`, and `src/shake_prng.c` →
`seedexpander_init`/`seedexpander`). Reproduce `y` in Python from
`SHAKE256(seed_dk ‖ 0x01)` using the **exact** C reduction and read pattern.
Once `y` matches, `x` will match (same call, next slice), and the rest of
keygen (`s = x + h·y`) can be checked against the intermediates' `s`.

A standalone Python reproduction script scaffold is in §6 to make this fast.

---

## 5. Summary of required code changes (once §4 is resolved)

| # | Change | File | Confidence |
|---|--------|------|-----------|
| 1 | Fix fixed-weight sampler to match reference exactly | `src/poly/sampling.rs` | **blocking, unknown** |
| 2 | `sigma` = K bytes, not 32 | `src/kem.rs` | high |
| 3 | Add KAT-format `sk` serializer (full 2321-B layout) | `src/parsing.rs` / `src/kem.rs` | high (layout TBD) |
| 4 | Add `aes` dev-dependency + NIST CTR_DRBG in harness | `Cargo.toml`, `tests/kat.rs` | high (DRBG verified) |
| 5 | Confirm DRBG draw order in keygen/encaps | `tests/kat.rs` | medium |
| 6 | Add `kat` feature, `.rsp` parser, byte-for-byte asserts | `Cargo.toml`, `tests/kat.rs` | high |

Also re-examine whether `x`/`y` should use the **same** sampler as `r1/r2/e`.
The current code uses rejection for `x/y` and Barrett for `r1/r2/e`; the
reference may use **one** sampler for all five (this is the long-standing
unverified assumption flagged in `CLAUDE.md` Step 14). Resolving §4 settles it.

---

## 6. Python scaffold to reproduce `y` (for the human to iterate)

```python
import hashlib, re
f = "tests/kat-vectors/hqc-128/intermediates_values"
text = open(f).read()
keygen = text.split('### KEYGEN ###')[1].split('### ENCAPS ###')[0]
def grab(l): 
    return bytes.fromhex(re.search(r'^'+l+r':\s*([0-9a-fA-F]+)\s*$', keygen, re.M).group(1))
seed_dk, y = grab('seed_dk'), grab('y')
N, OMEGA = 17669, 66
stream = hashlib.shake_256(seed_dk + b'\x01').digest(4000)   # verified XOF
yset = {i for i in range(N) if (y[i//8] >> (i % 8)) & 1}

# TODO: implement the EXACT reference vect_set_random_fixed_weight here,
# reading from `stream`, and assert set(positions) == yset.
```

Self-test for the DRBG side (pure-Python NIST CTR_DRBG, AES-256) lived in
`/tmp/aes_drbg.py` during the session; it validated against the canonical
`[0..47]` first seed. It should be ported into `tests/kat.rs` using the `aes`
crate once the sampler is fixed.

---

## 7. Bottom line

- `hash.rs`, `sample_uniform`, and the DRBG algorithm are **confirmed correct**.
- The **fixed-weight sampler is the blocker** — it does not match the reference,
  and the exact reference algorithm needs to be read from the 2025 C source.
- Secondary fixes (sigma length, full `sk` format, DRBG draw order) are
  understood and low-risk once the sampler is right.
- No `.rsp` byte-for-byte comparison can pass until §4 is fixed.
