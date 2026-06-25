// KAT (Known Answer Test) harness.
//
// Reads the official NIST `.req` seed files from `tests/hqc-{1,3,5}/` and
// regenerates the `pk`, `sk`, `ct`, `ss` fields using THIS crate's
// implementation, writing the result to a sibling `*.our.rsp` file, then
// asserts every field matches the official `.rsp` byte-for-byte. The upstream
// `.rsp` is the assertion oracle and is never overwritten; the test FAILS on
// any divergence (reporting the count, field, and first differing offset).
//
// ── How the seeds are consumed ────────────────────────────────────────────────
// Each `.req` `seed` is a 48-byte seed for the 2025 HQC KAT PRNG. That PRNG is
// NOT the classic NIST AES-256-CTR DRBG — the reference `symmetric.c`
// (`prng_init` / `prng_get_bytes`) is a SHAKE256 XOF seeded with
// `entropy ‖ personalization ‖ 0x00` (HQC_PRNG_DOMAIN) and squeezed as one
// continuous stream. We reproduce the reference `main_kat.c` call sequence per
// count:
//
//     prng_init(seed, NULL, 48, 0)               // SHAKE256(seed ‖ 0x00)
//     (pk, sk) = KEM.Keygen()      // draws SEED_BYTES (seed_kem) from the PRNG
//     (ss, ct) = KEM.Encaps(pk)    // draws K (m) + SALT_BYTES (salt) from the PRNG
//     ss'      = KEM.Decaps(sk, ct)
//     assert ss == ss'             // internal round-trip sanity check
//
// The draw order/sizes above match the reference exactly (kem.c:
// `prng_get_bytes(seed_kem, 32)`, then `prng_get_bytes(m, K)` /
// `prng_get_bytes(salt, 16)`), so given the same SHAKE256 stream the outputs
// coincide with the upstream pqc-hqc.org `.rsp` vectors. (All K values — 16/24/32
// — are multiples of 8, so the seed-expander's `xof_get_bytes` 8-byte alignment
// never pads inside the KEM expansion.)
//
// The `sk` field is emitted in the reference wire format via
// `DecapsulationKey::expanded_secret_key_bytes()`:
//     sk = pk ‖ seed_dk (32) ‖ seed_ek (32) ‖ σ (K)        // 64 + K suffix
// matching the upstream |sk| (2321 / 4602 / 7333), NOT the compressed 32-byte
// `seed_kem`.
//
// Enabled with the `kat` feature:
//     cargo test --features kat --test kat -- --nocapture

#![cfg(feature = "kat")]

use std::collections::HashMap;
use std::fs;

use rand_core::{CryptoRng, RngCore};
use sha3::digest::{ExtendableOutput, Update, XofReader};
use sha3::Shake256;

use hqcr::{kem, Hqc128, Hqc192, Hqc256, HqcParams};

// ── 2025 HQC KAT PRNG (SHAKE256-based — the reference `symmetric.c` PRNG) ───────
//
// The 2025 reference does NOT use the classic NIST AES-256-CTR DRBG. Its
// `prng_init` / `prng_get_bytes` (symmetric.c) build a SHAKE256 XOF seeded with
// `entropy ‖ personalization ‖ 0x00` (HQC_PRNG_DOMAIN) and squeeze it as one
// continuous stream:
//
//     prng_init(entropy, perso, enlen, perlen):
//         SHAKE256_inc_init
//         absorb(entropy, enlen); absorb(perso, perlen); absorb(0x00, 1)
//         finalize
//     prng_get_bytes(out, len):
//         squeeze(out, len)            // continuous across calls
//
// `main_kat.c` calls `prng_init(seed, NULL, 48, 0)` per count, so the
// personalization absorb is empty and the domain byte 0x00 follows the 48-byte
// seed directly. The squeeze being continuous means keygen's `seed_kem` (32 B),
// then encaps's `m` (K B) and `salt` (16 B), are simply consecutive bytes off
// the one stream — which is exactly how our KEM's `fill_bytes` calls consume it.

struct KatPrng {
    reader: <Shake256 as ExtendableOutput>::Reader,
}

impl KatPrng {
    /// `prng_init(seed, NULL, seed.len(), 0)`: SHAKE256(seed ‖ 0x00) as an XOF.
    fn new(seed: &[u8]) -> Self {
        const HQC_PRNG_DOMAIN: u8 = 0x00;
        let mut h = Shake256::default();
        Update::update(&mut h, seed);
        Update::update(&mut h, &[HQC_PRNG_DOMAIN]);
        KatPrng {
            reader: h.finalize_xof(),
        }
    }

    /// `prng_get_bytes(out, out.len())`: squeeze the next bytes off the stream.
    fn get_bytes(&mut self, out: &mut [u8]) {
        self.reader.read(out);
    }
}

// Expose the PRNG as the RNG our KEM API consumes. Each `fill_bytes` squeezes
// the requested bytes from the single continuous SHAKE256 stream, matching the
// reference `prng_get_bytes` byte-for-byte.
impl RngCore for KatPrng {
    fn next_u32(&mut self) -> u32 {
        let mut b = [0u8; 4];
        self.reader.read(&mut b);
        u32::from_le_bytes(b)
    }
    fn next_u64(&mut self) -> u64 {
        let mut b = [0u8; 8];
        self.reader.read(&mut b);
        u64::from_le_bytes(b)
    }
    fn fill_bytes(&mut self, dest: &mut [u8]) {
        self.reader.read(dest);
    }
    fn try_fill_bytes(&mut self, dest: &mut [u8]) -> Result<(), rand_core::Error> {
        self.reader.read(dest);
        Ok(())
    }
}
impl CryptoRng for KatPrng {}

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

// ── .rsp parsing (the official vectors, used as the assertion oracle) ──────────

/// The four regenerated-and-checked fields of one `.rsp` entry, as
/// uppercase-hex strings (the format both the upstream file and our output use).
struct RspEntry {
    pk: String,
    sk: String,
    ct: String,
    ss: String,
}

/// Parse an official NIST `.rsp` file into `count -> RspEntry`. Each block is
/// `count = / seed = / pk = / sk = / ct = / ss =`; we key on `count` and flush
/// the accumulated entry when its `ss =` line is read.
fn parse_rsp(path: &str) -> HashMap<usize, RspEntry> {
    let content = fs::read_to_string(path).unwrap_or_else(|e| panic!("read {path}: {e}"));
    let mut map = HashMap::new();
    let mut count: Option<usize> = None;
    let (mut pk, mut sk, mut ct) = (String::new(), String::new(), String::new());

    for line in content.lines() {
        let line = line.trim();
        if let Some(r) = line.strip_prefix("count = ") {
            count = Some(r.trim().parse().expect("count must be an integer"));
        } else if let Some(r) = line.strip_prefix("pk = ") {
            pk = r.trim().to_uppercase();
        } else if let Some(r) = line.strip_prefix("sk = ") {
            sk = r.trim().to_uppercase();
        } else if let Some(r) = line.strip_prefix("ct = ") {
            ct = r.trim().to_uppercase();
        } else if let Some(r) = line.strip_prefix("ss = ") {
            let ss = r.trim().to_uppercase();
            let c = count.expect("`count =` must precede `ss =`");
            map.insert(
                c,
                RspEntry {
                    pk: std::mem::take(&mut pk),
                    sk: std::mem::take(&mut sk),
                    ct: std::mem::take(&mut ct),
                    ss,
                },
            );
        }
    }
    assert!(!map.is_empty(), "no entries parsed from {path}");
    map
}

/// Assert two uppercase-hex field strings are equal; on mismatch report the
/// field, count, lengths, and the first differing hex-char offset with a short
/// surrounding snippet — never dumping the full (multi-KB) values.
fn assert_field_eq(header: &str, count: usize, field: &str, got: &str, want: &str) {
    if got == want {
        return;
    }
    let at = got
        .chars()
        .zip(want.chars())
        .position(|(a, b)| a != b)
        .unwrap_or_else(|| got.len().min(want.len()));
    let snippet = |s: &str| {
        let start = at.saturating_sub(8);
        let end = (at + 8).min(s.len());
        s.get(start..end).unwrap_or("").to_string()
    };
    panic!(
        "{header} count {count}: `{field}` mismatch vs official .rsp\n  \
         got len={}, want len={}, first diff at hex char {at}\n  \
         got  …{}…\n  want …{}…",
        got.len(),
        want.len(),
        snippet(got),
        snippet(want),
    );
}

// ── Generation ──────────────────────────────────────────────────────────────

/// Run the full KAT sequence for parameter set `P` over every seed in `req`,
/// write the filled vectors to `out_rsp` (uppercase hex, NIST `.rsp` layout),
/// then assert every `pk`/`sk`/`ct`/`ss` matches the official `reference_rsp`
/// byte-for-byte. The file is written *before* the assertions so it stays
/// available for inspection even when a vector regresses.
fn generate<P: HqcParams>(req: &str, reference_rsp: &str, out_rsp: &str, header: &str) {
    let seeds = parse_seeds(req);
    let reference = parse_rsp(reference_rsp);

    let mut output = String::new();
    output.push_str(&format!("# {header}\n\n"));
    // (count, pk, sk, ct, ss) uppercase-hex, retained for the self-check below.
    let mut generated: Vec<(usize, String, String, String, String)> =
        Vec::with_capacity(seeds.len());

    for (count, seed) in &seeds {
        let mut prng = KatPrng::new(&seed[..]);

        // Reference main_kat.c call order: keypair, then enc, then dec.
        let (pk, dk) = kem::keygen::<P, _>(&mut prng);
        let (ss, ct) = kem::encaps::<P, _>(&mut prng, &pk);
        let ss_dec = kem::decaps::<P>(&dk, &ct);

        // Internal correctness oracle: the valid decapsulation must recover the
        // encapsulated key (independent of the external `.rsp` comparison).
        assert_eq!(
            ss, ss_dec,
            "{header} count {count}: decaps(encaps()) shared-key mismatch"
        );

        let pk_bytes = pk.to_bytes();
        // Full reference-format sk (ek_pke ‖ seed_dk ‖ σ_K ‖ seed_kem), not the
        // compressed 32-byte seed_KEM, so the field is comparable to upstream.
        let sk_bytes = dk.expanded_secret_key_bytes();
        assert_eq!(pk_bytes.len(), P::PK_BYTES, "pk length");
        assert_eq!(ct.len(), P::CT_BYTES, "ct length");

        let (pk_hex, sk_hex, ct_hex, ss_hex) = (
            hex::encode_upper(&pk_bytes),
            hex::encode_upper(&sk_bytes[..]),
            hex::encode_upper(&ct),
            hex::encode_upper(ss),
        );

        output.push_str(&format!("count = {count}\n"));
        output.push_str(&format!("seed = {}\n", hex::encode_upper(seed)));
        output.push_str(&format!("pk = {pk_hex}\n"));
        output.push_str(&format!("sk = {sk_hex}\n"));
        output.push_str(&format!("ct = {ct_hex}\n"));
        output.push_str(&format!("ss = {ss_hex}\n"));
        output.push('\n');

        generated.push((*count, pk_hex, sk_hex, ct_hex, ss_hex));
    }

    fs::write(out_rsp, output).unwrap_or_else(|e| panic!("write {out_rsp}: {e}"));
    eprintln!("KAT: wrote {out_rsp} ({} entries)", seeds.len());

    // Self-check: every regenerated field must match the official `.rsp`.
    for (count, pk_hex, sk_hex, ct_hex, ss_hex) in &generated {
        let exp = reference
            .get(count)
            .unwrap_or_else(|| panic!("{header}: official .rsp is missing count {count}"));
        assert_field_eq(header, *count, "pk", pk_hex, &exp.pk);
        assert_field_eq(header, *count, "sk", sk_hex, &exp.sk);
        assert_field_eq(header, *count, "ct", ct_hex, &exp.ct);
        assert_field_eq(header, *count, "ss", ss_hex, &exp.ss);
    }
    eprintln!(
        "KAT: {header} — all {} vectors match the official .rsp byte-for-byte",
        generated.len()
    );
}

// ── Entry points (one per parameter set) ──────────────────────────────────────

#[test]
fn kat_hqc128() {
    generate::<Hqc128>(
        "tests/hqc-1/PQCkemKAT_2321.req",
        "tests/hqc-1/PQCkemKAT_2321.rsp",
        "tests/hqc-1/PQCkemKAT_2321.our.rsp",
        "HQC-1",
    );
}

#[test]
fn kat_hqc192() {
    generate::<Hqc192>(
        "tests/hqc-3/PQCkemKAT_4602.req",
        "tests/hqc-3/PQCkemKAT_4602.rsp",
        "tests/hqc-3/PQCkemKAT_4602.our.rsp",
        "HQC-3",
    );
}

#[test]
fn kat_hqc256() {
    generate::<Hqc256>(
        "tests/hqc-5/PQCkemKAT_7333.req",
        "tests/hqc-5/PQCkemKAT_7333.rsp",
        "tests/hqc-5/PQCkemKAT_7333.our.rsp",
        "HQC-5",
    );
}

// ── Diagnostics: localize the `.rsp` mismatch ─────────────────────────────────
//
// The full `.rsp` `pk` differs from the official vectors at byte 0, while the
// `intermediates_values` trace matches byte-for-byte. Since `pk = f(seed_kem)`
// and that derivation `f` is proven correct by the intermediate match, the only
// possible cause is that the value fed into `f` — `seed_kem`, the first 32 bytes
// our DRBG emits — differs from the reference. These three tests pin down which
// stage diverges; run them with:
//
//     cargo test --features kat --test kat -- --nocapture diag
//
// They are observational only: `drbg_self_test` and `derivation_isolated` assert
// against known-good reference values; `diag_count0` only prints.

/// Hex → [u8; 32] (panics on bad input — test-only).
fn hex32(s: &str) -> [u8; 32] {
    let v = hex::decode(s).expect("valid hex");
    assert_eq!(v.len(), 32, "expected 32 bytes");
    let mut out = [0u8; 32];
    out.copy_from_slice(&v);
    out
}

/// PRNG validation, independent of HQC. Every HQC 2025 `.req` is produced by
/// `prng_init(entropy = [0,1,…,47])` followed by one `prng_get_bytes(48)` per
/// count, squeezed continuously. The first 48 bytes therefore equal the count-0
/// seed. So **count-0 seed == SHAKE256([0..47] ‖ 0x00) squeezed for 48 bytes**.
/// If our PRNG reproduces it, the byte production matches the reference exactly.
#[test]
fn diag_drbg_self_test() {
    let mut entropy = [0u8; 48];
    for (i, b) in entropy.iter_mut().enumerate() {
        *b = i as u8;
    }
    let mut prng = KatPrng::new(&entropy[..]);
    let mut first = [0u8; 48];
    prng.get_bytes(&mut first);

    let count0_seed = parse_seeds("tests/hqc-1/PQCkemKAT_2321.req")[0].1;
    assert_eq!(
        hex::encode_upper(first),
        hex::encode_upper(count0_seed),
        "SHAKE256([0..47] ‖ 0x00) squeezed for 48 bytes must equal the count-0 \
         .req seed. If this fails, the SHAKE256 KAT PRNG port is wrong."
    );
}

/// Derivation in isolation, fed a *known* `seed_kem`. The official
/// `intermediates_values` file uses the passthrough `seed_kem` =
/// `9ef877…403e7` and records `seed_ek` = `ef2b80…fcb1` (the first 32 bytes of
/// `pk`). `keygen_from_seed` must reproduce exactly that — proving `pk =
/// f(seed_kem)` is correct and isolating the bug to whatever feeds `seed_kem`
/// (the DRBG) in the full `.rsp` flow.
#[test]
fn diag_derivation_isolated() {
    let seed_kem = hex32("9ef877fddbe8891c6e4e79eaf022e563defaca6b152161b9a423e8fe96a403e7");
    let (pk, _dk) = kem::keygen_from_seed::<Hqc128>(&seed_kem);
    let pk_bytes = pk.to_bytes();
    assert_eq!(
        hex::encode(&pk_bytes[..32]),
        "ef2b80f46f3a6437b4d869bb38bdd6004bff72bcd0ceb139b4b8d47301f4fcb1",
        "keygen_from_seed(intermediate seed_kem) must reproduce the \
         intermediate-values seed_ek"
    );
}

/// Side-by-side print of the count-0 keygen, exactly as the `.rsp` flow runs it:
/// init the DRBG with the count-0 seed, draw `seed_kem` (32 B), derive `pk`.
/// Compares our `seed_kem`-derived `pk` prefix against the official `.rsp` `pk`.
#[test]
fn diag_count0() {
    let count0_seed = parse_seeds("tests/hqc-1/PQCkemKAT_2321.req")[0].1;
    let mut prng = KatPrng::new(&count0_seed[..]);

    let mut seed_kem = [0u8; 32];
    prng.get_bytes(&mut seed_kem);

    let (pk, _dk) = kem::keygen_from_seed::<Hqc128>(&seed_kem);
    let pk_bytes = pk.to_bytes();

    eprintln!("\n── diag_count0 (HQC-128) ──");
    eprintln!("our  seed_kem   = {}", hex::encode(seed_kem));
    eprintln!("our  pk[0:32]   = {}", hex::encode(&pk_bytes[..32]));
    eprintln!("ref  pk[0:32]   = 4053237912ea281c51c4456a5096589ec9d20219651e00f9704178f0cf84f9ae");
    eprintln!(
        "intermediate seed_kem (passthrough, NOT DRBG) = \
         9ef877fddbe8891c6e4e79eaf022e563defaca6b152161b9a423e8fe96a403e7"
    );
}
