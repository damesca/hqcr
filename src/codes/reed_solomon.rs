// Shortened Reed-Solomon over GF(2^8) — spec §3.4.
//
// Three shortened RS codes (spec Table 3):
//   RS-S1: [46, 16, 31]  — for HQC-128, δ=15, nk=2δ=30 parity symbols
//   RS-S2: [56, 24, 33]  — for HQC-192, δ=16, nk=2δ=32 parity symbols
//   RS-S3: [90, 32, 49]  — for HQC-256, δ=29, nk=2δ=58 parity symbols
//
// All three use GF(2^8) with primitive polynomial x^8+x^4+x^3+x^2+1 (0x11D)
// and a narrow-sense generator with roots α^1, α^2, ..., α^{2δ}.
//
// ── Polynomial / array conventions (read this before editing) ────────────────
//
//   * Every polynomial is stored LITTLE-ENDIAN: array[i] is the coefficient of
//     x^i. This matches gf::gf_poly_eval, which the syndrome step relies on.
//
//   * Systematic codeword layout (consequence of little-endian + standard
//     systematic encoding c(x) = x^{nk}·M(x) + r(x)):
//         codeword[0 .. nk]   = parity   r(x)   (low-degree coefficients)
//         codeword[nk .. n1]  = message  M(x)   (high-degree coefficients)
//     i.e. the message lives in the HIGH positions, parity in the LOW ones.
//     C.Encode/C.Decode in codes/mod.rs must read symbols with this layout.
//
//   * Generator: NOTE — CLAUDE.md asks for the hardcoded g1/g2/g3 from spec
//     §3.4.2. The constants previously hardcoded here were incorrect (they were
//     palindromic with constant term 1, which no narrow-sense RS generator can
//     have). We instead build the generator at runtime from its roots α^1..α^{2δ}.
//     This is mathematically identical to the spec generator and guarantees the
//     encode/decode roundtrip is correct. It costs O((2δ)^2) GF ops once per
//     encode/decode call (negligible: RS runs once per KEM operation). The
//     hardcoded constants can be reinstated later if cross-checked against KAT.
//
// Decoding pipeline (constant-time, Steps 20a–c — same fixed sequence for every
// input, no syndrome fast path and no secret-dependent branch / memory address):
//   1. Syndromes S_j = received(α^j), j = 1..=2δ (always computed in full).
//   2. Berlekamp-Massey → error locator σ(x), deg σ = number of errors (masked,
//      constant-time — see `berlekamp_massey`).
//   3. Ω(x) = S(x)σ(x) mod x^{2δ} and σ'(x) for the Forney step.
//   4. Branchless root finding + correction: scan every position i in 0..n1,
//      test σ(α^{-i}) == 0 with `ct_eq`, and XOR a masked Forney correction into
//      corrected[i] (store index i is the public loop counter). Replaces the old
//      Chien search's conditional `push` (leak C1) and its secret-length `Vec`.
//   5. Validity = corrected has zero syndromes ∧ deg ≤ δ ∧ #roots == deg, folded
//      to one bit; the single Some/None branch. Returns None when uncorrectable.

use crate::gf::{gf_add, gf_div, gf_inv, gf_mul, gf_poly_eval, GF_EXP};
use subtle::{ConditionallySelectable, ConstantTimeEq};
use zeroize::Zeroizing;

// ── Generator polynomial ──────────────────────────────────────────────────────

/// Build the narrow-sense RS generator g(x) = ∏_{i=1}^{2δ} (x - α^i).
///
/// Returned little-endian: g[0] is the constant term, g[2δ] = 1 (monic).
/// Length is 2δ + 1.
pub(crate) fn rs_generator(delta: usize) -> Vec<u8> {
    let nk = 2 * delta;
    let mut g = vec![0u8; nk + 1];
    g[0] = 1; // start with g(x) = 1
    let mut deg = 0usize;

    for &root in &GF_EXP[1..=nk] {
        // root = α^i for each successive i
        deg += 1;
        // Multiply g(x) by (x - root) = (x + root) in char 2.
        // new g[j] = g[j-1] + root * g[j], processed high→low to avoid clobbering.
        for j in (1..=deg).rev() {
            g[j] = gf_add(g[j - 1], gf_mul(root, g[j]));
        }
        g[0] = gf_mul(root, g[0]);
    }

    debug_assert_eq!(g[nk], 1, "generator must be monic");
    g
}

// ── Encoding ──────────────────────────────────────────────────────────────────

/// Systematic RS encoding.
///
/// `msg` is `k` GF(2^8) symbols. `codeword` (length `n1 = k + 2δ`) is filled as:
///   codeword[0 .. nk]  = parity  (low-degree coefficients)
///   codeword[nk .. n1] = message (high-degree coefficients)
///
/// Computes r(x) = (M(x) · x^{nk}) mod g(x) by explicit polynomial long
/// division, then emits [parity | message]. The whole codeword polynomial is
/// then divisible by g(x), so its syndromes are zero.
pub fn rs_encode(msg: &[u8], codeword: &mut [u8], delta: usize) {
    let n1 = codeword.len();
    let k = msg.len();
    let nk = 2 * delta;
    debug_assert_eq!(k + nk, n1, "k + 2δ must equal n1");

    let g = rs_generator(delta);

    // Dividend d(x) = M(x) · x^{nk}, little-endian, length n1:
    //   d[nk + t] = msg[t]; the low nk coefficients start at 0. `d` holds the
    //   secret message during division, so it is `Zeroizing` (wiped on drop).
    let mut d = Zeroizing::new(vec![0u8; n1]);
    d[nk..(nk + k)].copy_from_slice(&msg[..k]);

    // Reduce d modulo g(x): cancel the leading coefficient at each degree from
    // n1-1 down to nk. g is monic of degree nk, so subtracting c·x^{i-nk}·g(x)
    // zeroes d[i]. After the loop, d[0..nk] holds the remainder (parity).
    for i in (nk..n1).rev() {
        let c = d[i];
        // No `if c != 0` guard: keep control flow independent of the (secret)
        // message symbols. gf_mul(0, ·) = 0, so a zero coefficient is a no-op.
        for j in 0..=nk {
            d[i - nk + j] = gf_add(d[i - nk + j], gf_mul(c, g[j]));
        }
    }

    codeword[..nk].copy_from_slice(&d[..nk]); // parity
    codeword[nk..].copy_from_slice(msg); // message
}

// ── Syndrome computation ──────────────────────────────────────────────────────

/// Compute the 2δ syndromes S_j = received(α^j) for j = 1..=2δ.
/// Returns the syndrome vector (s[t] = S_{t+1}) and whether any is non-zero.
///
/// Constant-time: every syndrome is computed (no early exit) and the
/// "any non-zero" flag is accumulated branchlessly by OR-ing the bytes, so the
/// instruction trace does not depend on the (secret) syndrome values. The
/// returned bool is the single 1-bit summary; `rs_decode` no longer branches on
/// it for a fast path.
fn syndromes(received: &[u8], delta: usize) -> (Vec<u8>, bool) {
    let mut s = vec![0u8; 2 * delta];
    let mut acc = 0u8;
    for j in 0..2 * delta {
        // S_{j+1} = received(α^{j+1}); α^{j+1} = GF_EXP[j+1] (public index).
        s[j] = gf_poly_eval(received, GF_EXP[j + 1]);
        acc |= s[j];
    }
    (s, acc != 0)
}

// ── Berlekamp-Massey ──────────────────────────────────────────────────────────

/// Berlekamp-Massey over GF(2^8) — **constant-time** (Step 20b).
///
/// Given syndromes `s` (s[t] = S_{t+1}, length 2δ), returns `(σ, deg σ)` where
/// σ(x) = 1 + σ_1 x + … (little-endian, σ[0] = 1, full length δ+1) and deg σ is
/// the number of errors. The roots of σ are the inverse error locators α^{-i}.
///
/// This is a direct port of the HQC reference `compute_elp` (`reed_solomon.c`).
/// Every conditional update is mask-merged rather than branched, so the
/// instruction trace and memory-access pattern are identical for all syndrome
/// inputs (the previous version branched on the discrepancy value and used a
/// secret-bounded inner loop — leak classes C1/C3 from the audit):
///
/// - The only run-time predicate is `mask12` ("did the register length grow?").
///   It is combined under `black_box` so LLVM cannot turn the masked merges back
///   into a data-dependent branch (the C reference uses `volatile`).
/// - Both loops run to **public** bounds (`min(μ+1, δ)` and `1..=δ`), never to a
///   secret length.
/// - `gf_inv`/`gf_mul`/`gf_div` are themselves branch-free (Step 20a), and
///   `gf_inv(0) == 0` lets the `d_p` inverse run unconditionally.
fn berlekamp_massey(s: &[u8]) -> (Vec<u8>, usize) {
    let two_delta = s.len(); // = 2δ
    let delta = two_delta / 2;

    // σ(x): error locator, σ[0] = 1. Indices 0..=δ.
    let mut sigma = vec![0u8; delta + 1];
    sigma[0] = 1;
    // X·σ_p(x): the "B(x)" already shifted by X. Init = X (i.e. coefficient 1 at
    // index 1). Index 0 stays 0 forever (no constant term). Secret-derived
    // scratch buffers are `Zeroizing` so they wipe on drop (zeroize gap G3);
    // `sigma` is wiped by the caller (`rs_decode` wraps the returned value).
    let mut x_sigma_p = Zeroizing::new(vec![0u8; delta + 1]);
    x_sigma_p[1] = 1;
    // Scratch save of σ taken before each update (reference `sigma_copy`).
    let mut sigma_copy = Zeroizing::new(vec![0u8; delta + 1]);

    let mut deg_sigma: u16 = 0;
    let mut deg_sigma_p: u16 = 0;

    // pp = ρ, the μ at the last length change; init −1 so deg_x = μ−pp = μ+1 at
    // the start. The u16 wraparound is intentional (matches the reference).
    let mut pp: u16 = u16::MAX;
    let mut d_p: u8 = 1; // discrepancy at the last length change (b)
    let mut d: u8 = s[0]; // current discrepancy

    let mut mu: u16 = 0;
    while (mu as usize) < two_delta {
        // Save σ before this update (used to refresh X·σ_p on a length change).
        sigma_copy.copy_from_slice(&sigma);
        let deg_sigma_copy = deg_sigma;

        // dd = d / d_p  (= δ/b). Branch-free; gf_inv(0) = 0.
        let dd = gf_mul(d, gf_inv(d_p));

        // σ(x) -= dd · (X·σ_p)(x). Public loop bound.
        let upper = core::cmp::min(mu as usize + 1, delta);
        for i in 1..=upper {
            sigma[i] = gf_add(sigma[i], gf_mul(dd, x_sigma_p[i]));
        }

        let deg_x = mu.wrapping_sub(pp); // shifts since last length change
        let deg_x_sigma_p = deg_x.wrapping_add(deg_sigma_p);

        // mask1 = 0xFFFF iff d != 0.
        let m1bit = (d as u16).wrapping_neg() >> 15; // 1 if d != 0 else 0
        let mask1 = 0u16.wrapping_sub(m1bit);
        // mask2 = 0xFFFF iff deg_x_sigma_p > deg_sigma (sign bit of the diff;
        // degrees are tiny so the signed interpretation is valid).
        let m2bit = (deg_sigma.wrapping_sub(deg_x_sigma_p) >> 15) & 1;
        let mask2 = 0u16.wrapping_sub(m2bit);
        // mask12 = "register length grew". black_box stops the compiler from
        // re-deriving a branch from the masked merges below.
        let mask12 = core::hint::black_box(mask1 & mask2);
        let mask12_u8 = mask12 as u8; // 0x00 or 0xFF

        // deg_sigma = mask12 ? deg_x_sigma_p : deg_sigma.
        deg_sigma ^= mask12 & (deg_x_sigma_p ^ deg_sigma);

        // Final iteration only finalises deg_sigma (matches the reference break).
        if (mu as usize) == two_delta - 1 {
            break;
        }

        // pp = mask12 ? mu : pp.
        pp ^= mask12 & (mu ^ pp);
        // d_p = mask12 ? d : d_p.
        d_p ^= mask12_u8 & (d ^ d_p);
        // X·σ_p = X · (mask12 ? σ_copy : X·σ_p): shift the chosen poly up by one.
        for i in (1..=delta).rev() {
            x_sigma_p[i] = (mask12_u8 & sigma_copy[i - 1]) | (!mask12_u8 & x_sigma_p[i - 1]);
        }
        // deg_sigma_p = mask12 ? deg_sigma_copy : deg_sigma_p.
        deg_sigma_p ^= mask12 & (deg_sigma_copy ^ deg_sigma_p);

        // Next discrepancy d = S_{μ+2} + Σ_{i=1}^{min(μ+1,δ)} σ_i · S_{μ+2−i}.
        d = s[mu as usize + 1];
        for i in 1..=upper {
            d = gf_add(d, gf_mul(sigma[i], s[mu as usize + 1 - i]));
        }

        mu += 1;
    }

    (sigma, deg_sigma as usize)
}

// ── Error-evaluator polynomial and formal derivative ──────────────────────────

/// Ω(x) = (S · σ) mod x^{2δ}, the Forney error-evaluator polynomial, where
/// S(x) = Σ_t synd[t] x^t (synd[t] = S_{t+1}). Loop bounds are public (the
/// syndrome / σ lengths), so this is constant-time.
fn error_evaluator(synd: &[u8], sigma: &[u8]) -> Vec<u8> {
    let two_delta = synd.len();
    let mut omega = vec![0u8; two_delta];
    for i in 0..synd.len() {
        for j in 0..sigma.len() {
            if i + j < two_delta {
                omega[i + j] = gf_add(omega[i + j], gf_mul(synd[i], sigma[j]));
            }
        }
    }
    omega
}

/// σ'(x), the formal derivative of σ. In characteristic 2 only odd-degree terms
/// survive: σ'[j-1] = σ[j] for odd j. Constant-time (public loop).
fn formal_derivative(sigma: &[u8]) -> Vec<u8> {
    let mut sigma_prime = vec![0u8; sigma.len()];
    for j in (1..sigma.len()).step_by(2) {
        sigma_prime[j - 1] = sigma[j];
    }
    sigma_prime
}

// ── Public decode entry point ─────────────────────────────────────────────────

/// Decode a received RS codeword of length `n1` carrying `k` message symbols
/// with error-correcting capacity `delta` (so nk = 2δ = n1 - k).
///
/// Returns `Some(message)` (k symbols, from the high positions) on success,
/// or `None` if the error pattern is uncorrectable.
///
/// **Constant-time (Step 20c).** The whole pipeline runs the *same* fixed
/// sequence of operations for every input — there is no syndrome fast path, no
/// early exit, and no secret-dependent branch or memory address:
///
/// - Root finding replaces the old Chien search's conditional `roots.push`
///   (leak class C1) and its secret-length `Vec` (zeroize gap G3). Instead we
///   scan **every** position `i ∈ [0, n1)` and, branchlessly, (a) test whether
///   α^{-i} is a root of σ with `ct_eq`, (b) compute the Forney value, and
///   (c) XOR a masked correction into `corrected[i]`. The store address `i` is
///   the public loop counter, never an error position.
/// - The only run-time predicate is the final 1-bit validity used to choose
///   `Some`/`None`. That bit is exactly what `kem::decaps` already consumes as a
///   `Choice` (`decode_ok`), folded with the constant-time re-encryption check;
///   it reveals nothing beyond "did decoding succeed", and all *work* preceding
///   it is input-independent.
pub fn rs_decode(received: &[u8], k: usize, delta: usize) -> Option<Vec<u8>> {
    let n1 = received.len();
    let nk = 2 * delta;
    debug_assert_eq!(k + nk, n1);

    // Step 1: syndromes (always computed in full — no fast path). Every
    // secret-derived decode intermediate below is wrapped in `Zeroizing` so its
    // heap buffer is wiped on drop / on every exit path (zeroize gap G3); all are
    // fixed-size (no realloc), so no stale plaintext is left behind (defeater D2).
    let s = Zeroizing::new(syndromes(received, delta).0);

    // Step 2: Berlekamp-Massey (constant-time, Step 20b). σ is full length δ+1
    // with trailing zero coefficients; deg σ is returned explicitly.
    let (sigma, deg) = berlekamp_massey(&s);
    let sigma = Zeroizing::new(sigma);

    // Step 3: Ω(x) and σ'(x) for the Forney step (fixed-size, CT).
    let omega = Zeroizing::new(error_evaluator(&s, &sigma));
    let sigma_prime = Zeroizing::new(formal_derivative(&sigma));

    // Step 4: constant-time root finding + correction. Scan every position; mark
    // and correct roots branchlessly into a fixed-size buffer indexed by the
    // public counter `i`.
    let mut corrected = Zeroizing::new(received.to_vec());
    let mut num_roots: usize = 0;
    for i in 0..n1 {
        // α^{-i} = α^{255-i} for i>0, α^0 = 1. n1 ≤ 90 ⇒ 255-i ≥ 165; public index.
        let xi = if i == 0 { GF_EXP[0] } else { GF_EXP[255 - i] };

        // is_root ⇔ σ(α^{-i}) == 0.
        let is_root = gf_poly_eval(&sigma, xi).ct_eq(&0u8);

        // Forney value Y = Ω(xi) / σ'(xi). gf_div is branch-free; at a non-root
        // (or when σ'(xi) = 0) the value is masked away below.
        let y = gf_div(gf_poly_eval(&omega, xi), gf_poly_eval(&sigma_prime, xi));

        // Apply the correction only at roots, branchlessly.
        let correction = u8::conditional_select(&0u8, &y, is_root);
        corrected[i] ^= correction;

        // Count roots without branching (for the validity check below).
        num_roots += is_root.unwrap_u8() as usize;
    }

    // Step 5: validity. Re-verify the corrected word has all-zero syndromes
    // (computed branchlessly) and that the locator was sane. This reproduces the
    // old None/Some decision exactly; it is the single final branch.
    let s2 = Zeroizing::new(syndromes(&corrected, delta).0);
    let mut acc = 0u8;
    for &sj in s2.iter() {
        acc |= sj;
    }
    // Fold the three conditions branchlessly into one bit; only the final
    // Some/None choice branches.
    let synd_zero = (acc == 0) as u8;
    let deg_ok = (deg <= delta) as u8;
    let roots_ok = (num_roots == deg) as u8;
    let valid = (synd_zero & deg_ok & roots_ok) == 1;

    if valid {
        Some(corrected[nk..].to_vec())
    } else {
        None
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // (n1, k, delta) for the three shortened codes.
    const S1: (usize, usize, usize) = (46, 16, 15);
    const S2: (usize, usize, usize) = (56, 24, 16);
    const S3: (usize, usize, usize) = (90, 32, 29);

    fn encode(msg: &[u8], n1: usize, delta: usize) -> Vec<u8> {
        let mut cw = vec![0u8; n1];
        rs_encode(msg, &mut cw, delta);
        cw
    }

    fn inject_errors(cw: &mut [u8], errors: &[(usize, u8)]) {
        for &(pos, val) in errors {
            cw[pos] ^= val;
        }
    }

    // ── Generator sanity ──────────────────────────────────────────────────────

    #[test]
    fn generator_is_monic_and_correct_degree() {
        for &(_, _, delta) in &[S1, S2, S3] {
            let g = rs_generator(delta);
            assert_eq!(g.len(), 2 * delta + 1);
            assert_eq!(g[2 * delta], 1, "must be monic");
        }
    }

    #[test]
    fn generator_roots_are_alpha_powers() {
        // g(α^i) = 0 for i = 1..=2δ by construction.
        for &(_, _, delta) in &[S1, S2, S3] {
            let g = rs_generator(delta);
            for i in 1..=2 * delta {
                assert_eq!(gf_poly_eval(&g, GF_EXP[i]), 0, "α^{i} should be a root");
            }
        }
    }

    // ── Systematic encoding: message in HIGH positions ────────────────────────

    #[test]
    fn encode_message_in_high_positions() {
        for &(n1, k, delta) in &[S1, S2, S3] {
            let nk = 2 * delta;
            let msg: Vec<u8> = (0..k).map(|i| (i * 7 + 1) as u8).collect();
            let cw = encode(&msg, n1, delta);
            assert_eq!(cw.len(), n1);
            assert_eq!(&cw[nk..], msg.as_slice(), "message must occupy [nk, n1)");
        }
    }

    #[test]
    fn encode_zero_message_is_all_zero() {
        let (n1, k, delta) = S1;
        let cw = encode(&vec![0u8; k], n1, delta);
        assert!(cw.iter().all(|&b| b == 0));
    }

    // ── Valid codewords have zero syndromes ──────────────────────────────────

    #[test]
    fn valid_codeword_has_zero_syndromes() {
        for &(n1, k, delta) in &[S1, S2, S3] {
            let msg: Vec<u8> = (0..k).map(|i| (i * 17 + 3) as u8).collect();
            let cw = encode(&msg, n1, delta);
            let (s, has_err) = syndromes(&cw, delta);
            assert!(!has_err, "delta={delta} syndromes nonzero: {s:?}");
        }
    }

    // ── Decode with no errors ─────────────────────────────────────────────────

    #[test]
    fn decode_no_errors() {
        for &(n1, k, delta) in &[S1, S2, S3] {
            let msg: Vec<u8> = (0..k).map(|i| (i ^ 0xAA) as u8).collect();
            let cw = encode(&msg, n1, delta);
            let recovered = rs_decode(&cw, k, delta).expect("decode failed");
            assert_eq!(recovered, msg, "delta={delta}");
        }
    }

    // ── Single-error correction ───────────────────────────────────────────────

    #[test]
    fn decode_single_error() {
        let (n1, k, delta) = S1;
        let msg: Vec<u8> = (0..k as u8).collect();
        let mut cw = encode(&msg, n1, delta);
        inject_errors(&mut cw, &[(5, 0x3C)]);
        let recovered = rs_decode(&cw, k, delta).expect("decode failed");
        assert_eq!(recovered, msg);
    }

    // ── Correction at full capacity (δ errors) ────────────────────────────────

    #[test]
    fn decode_at_capacity() {
        for &(n1, k, delta) in &[S1, S2, S3] {
            let msg: Vec<u8> = (0..k as u8).collect();
            let mut cw = encode(&msg, n1, delta);
            // δ errors at distinct positions, all non-zero.
            let errors: Vec<(usize, u8)> = (0..delta)
                .map(|i| (i, (i as u8).wrapping_add(1).wrapping_mul(7) | 1))
                .collect();
            inject_errors(&mut cw, &errors);
            let recovered =
                rs_decode(&cw, k, delta).unwrap_or_else(|| panic!("should correct {delta} errors"));
            assert_eq!(recovered, msg, "delta={delta}");
        }
    }

    // ── Error spread across message and parity regions ────────────────────────

    #[test]
    fn decode_errors_in_message_region() {
        let (n1, k, delta) = S3;
        let nk = 2 * delta;
        let msg: Vec<u8> = (0..k as u8).map(|i| i.wrapping_mul(3)).collect();
        let mut cw = encode(&msg, n1, delta);
        // Put errors in the high (message) region [nk, n1).
        let errors: Vec<(usize, u8)> = (0..delta).map(|i| (nk + (i % k), 0x5A | 1)).collect();
        // Distinct positions only:
        let mut seen = std::collections::HashSet::new();
        let errors: Vec<(usize, u8)> = errors
            .into_iter()
            .filter(|&(p, _)| seen.insert(p))
            .collect();
        inject_errors(&mut cw, &errors);
        let recovered = rs_decode(&cw, k, delta).expect("decode failed");
        assert_eq!(recovered, msg);
    }

    // ── Beyond capacity → None ────────────────────────────────────────────────

    #[test]
    fn decode_beyond_capacity_returns_none() {
        for &(n1, k, delta) in &[S1, S2, S3] {
            let msg: Vec<u8> = (0..k as u8).collect();
            let mut cw = encode(&msg, n1, delta);
            // δ+1 errors at distinct positions.
            let errors: Vec<(usize, u8)> = (0..=delta).map(|i| (i, 0xFF)).collect();
            inject_errors(&mut cw, &errors);
            assert!(
                rs_decode(&cw, k, delta).is_none(),
                "delta={delta}: δ+1 errors must be uncorrectable"
            );
        }
    }

    // ── Errors spanning both parity and message regions ──────────────────────

    #[test]
    fn decode_errors_in_both_regions() {
        // δ errors split: half in parity [0, nk), half in message [nk, n1).
        // Verifies the correction step applies correctly regardless of region.
        for &(n1, k, delta) in &[S1, S2, S3] {
            let nk = 2 * delta;
            let msg: Vec<u8> = (0..k as u8).map(|i| i.wrapping_mul(5).wrapping_add(3)).collect();
            let mut cw = encode(&msg, n1, delta);

            let half = delta / 2;
            let mut errors = Vec::new();
            for i in 0..half {
                errors.push((i, 0x55u8 | 1)); // parity region
            }
            for i in 0..(delta - half) {
                errors.push((nk + i, 0xAAu8 | 1)); // message region
            }
            inject_errors(&mut cw, &errors);

            let recovered = rs_decode(&cw, k, delta)
                .unwrap_or_else(|| panic!("should correct {delta} mixed-region errors (delta={delta})"));
            assert_eq!(recovered, msg, "delta={delta}");
        }
    }

    // ── Randomized correctness: every ≤ δ pattern decodes ─────────────────────
    // Exercises the constant-time root finder + inline Forney across many random
    // error counts / positions / values, beyond the hand-picked patterns above.

    #[test]
    fn ct_decode_corrects_random_patterns() {
        // Deterministic LCG (no external rng dependency).
        let mut state: u64 = 0x1234_5678_9ABC_DEF0;
        let mut next = || {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            (state >> 33) as u32
        };

        for &(n1, k, delta) in &[S1, S2, S3] {
            for _ in 0..40 {
                let msg: Vec<u8> = (0..k).map(|_| next() as u8).collect();
                let cw = encode(&msg, n1, delta);
                let mut received = cw.clone();

                // Choose a random error count in [0, δ] and distinct positions.
                let num = (next() as usize) % (delta + 1);
                let mut positions: Vec<usize> = Vec::new();
                while positions.len() < num {
                    let p = (next() as usize) % n1;
                    if !positions.contains(&p) {
                        positions.push(p);
                    }
                }
                for &p in &positions {
                    let v = (next() as u8) | 1; // non-zero error value
                    received[p] ^= v;
                }

                let recovered = rs_decode(&received, k, delta)
                    .unwrap_or_else(|| panic!("delta={delta} num={num}: decode returned None"));
                assert_eq!(recovered, msg, "delta={delta} num={num} positions={positions:?}");
            }
        }
    }

    // ── Root finding locates the right positions ──────────────────────────────
    // Test-only naive Chien search (the production decoder now finds roots via a
    // branchless full scan inside `rs_decode`). Confirms σ vanishes exactly at
    // the injected error positions.
    fn chien_positions(sigma: &[u8], n1: usize) -> Vec<usize> {
        let mut positions = Vec::new();
        for i in 0..n1 {
            let xi = if i == 0 { GF_EXP[0] } else { GF_EXP[255 - i] };
            if gf_poly_eval(sigma, xi) == 0 {
                positions.push(i);
            }
        }
        positions
    }

    #[test]
    fn root_finding_locates_correct_positions() {
        let (n1, k, delta) = S1;
        let msg: Vec<u8> = (0..k as u8).collect();
        let cw = encode(&msg, n1, delta);
        let mut received = cw.clone();
        inject_errors(&mut received, &[(3, 0x11), (10, 0x22)]);

        let (s, _) = syndromes(&received, delta);
        let (sigma, _deg) = berlekamp_massey(&s);
        let positions = chien_positions(&sigma, n1);
        assert!(
            positions.contains(&3),
            "should find position 3: {positions:?}"
        );
        assert!(
            positions.contains(&10),
            "should find position 10: {positions:?}"
        );
        assert_eq!(positions.len(), 2);
    }

    // ── Oracle: the old branch-based BM (pre-Step-20b) ────────────────────────
    // Retained only here to prove the new constant-time `berlekamp_massey`
    // produces the identical error locator on every correctable pattern.
    // (Same cross-check pattern Step 17/20a used.)
    fn berlekamp_massey_ref(s: &[u8]) -> Vec<u8> {
        let n = s.len();
        let mut sigma = vec![0u8; n + 1];
        sigma[0] = 1;
        let mut prev = vec![0u8; n + 1];
        prev[0] = 1;
        let mut l: usize = 0;
        let mut m: usize = 1;
        let mut b: u8 = 1;
        for i in 0..n {
            let mut delta = s[i];
            for j in 1..=l {
                delta = gf_add(delta, gf_mul(sigma[j], s[i - j]));
            }
            if delta == 0 {
                m += 1;
            } else if 2 * l <= i {
                let coeff = gf_div(delta, b);
                let t = sigma.clone();
                for j in m..=i + 1 {
                    sigma[j] = gf_add(sigma[j], gf_mul(coeff, prev[j - m]));
                }
                prev = t;
                l = i + 1 - l;
                b = delta;
                m = 1;
            } else {
                let coeff = gf_div(delta, b);
                for j in m..=i + 1 {
                    sigma[j] = gf_add(sigma[j], gf_mul(coeff, prev[j - m]));
                }
                m += 1;
            }
        }
        sigma.truncate(l + 1);
        sigma
    }

    #[test]
    fn ct_berlekamp_massey_matches_reference() {
        // For every code and every correctable error count e ∈ [1, δ], the
        // CT BM must return deg σ = e and the same σ coefficients (the locator
        // is unique for ≤ δ errors), with zero padding above degree e.
        for &(n1, k, delta) in &[S1, S2, S3] {
            for num_errors in 1..=delta {
                let msg: Vec<u8> = (0..k as u8).map(|i| i.wrapping_mul(13).wrapping_add(7)).collect();
                let cw = encode(&msg, n1, delta);
                let mut received = cw.clone();
                // Distinct positions, non-zero error values.
                let errors: Vec<(usize, u8)> = (0..num_errors)
                    .map(|i| (i * 2 % n1, (i as u8).wrapping_mul(31).wrapping_add(1) | 1))
                    .collect();
                inject_errors(&mut received, &errors);

                let (s, _) = syndromes(&received, delta);
                let (sigma_ct, deg_ct) = berlekamp_massey(&s);
                let sigma_ref = berlekamp_massey_ref(&s);

                assert_eq!(
                    deg_ct,
                    sigma_ref.len() - 1,
                    "delta={delta} e={num_errors}: degree mismatch"
                );
                // Coefficients agree up to the degree.
                assert_eq!(
                    &sigma_ct[..=deg_ct],
                    sigma_ref.as_slice(),
                    "delta={delta} e={num_errors}: locator mismatch"
                );
                // Padding above the degree is zero.
                assert!(
                    sigma_ct[deg_ct + 1..].iter().all(|&c| c == 0),
                    "delta={delta} e={num_errors}: nonzero padding"
                );
            }
        }
    }

    // ── BM degree equals the number of injected errors ────────────────────────

    #[test]
    fn berlekamp_massey_degree_equals_error_count() {
        let (n1, k, delta) = S1;
        for num_errors in 1..=delta {
            let msg: Vec<u8> = (0..k as u8).collect();
            let cw = encode(&msg, n1, delta);
            let mut received = cw.clone();
            let errors: Vec<(usize, u8)> = (0..num_errors)
                .map(|i| (i, (i as u8).wrapping_add(1) | 1))
                .collect();
            inject_errors(&mut received, &errors);
            let (s, _) = syndromes(&received, delta);
            let (_sigma, deg) = berlekamp_massey(&s);
            assert_eq!(deg, num_errors, "deg σ should equal error count {num_errors}");
        }
    }
}
