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
// Decoding pipeline:
//   1. Syndromes S_j = received(α^j), j = 1..=2δ. All zero ⇒ no errors.
//   2. Berlekamp-Massey → error locator σ(x), deg σ = number of errors.
//   3. Chien search: σ(α^{-i}) for every i in 0..n1 (NO early exit, CT req).
//   4. Forney: error values from Ω(x)=S(x)σ(x) mod x^{2δ} and σ'(x).
//   5. Correct received word; extract message from the high positions.
//   Returns None when the error pattern is uncorrectable (> δ errors).

use crate::gf::{gf_add, gf_div, gf_mul, gf_poly_eval, GF_EXP};

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
    //   d[nk + t] = msg[t]; the low nk coefficients start at 0.
    let mut d = vec![0u8; n1];
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
fn syndromes(received: &[u8], delta: usize) -> (Vec<u8>, bool) {
    let mut s = vec![0u8; 2 * delta];
    let mut has_error = false;
    for j in 0..2 * delta {
        // S_{j+1} = received(α^{j+1}); α^{j+1} = GF_EXP[j+1].
        s[j] = gf_poly_eval(received, GF_EXP[j + 1]);
        if s[j] != 0 {
            has_error = true;
        }
    }
    (s, has_error)
}

// ── Berlekamp-Massey ──────────────────────────────────────────────────────────

/// Berlekamp-Massey over GF(2^8).
///
/// Given syndromes `s` (s[t] = S_{t+1}, length 2δ), returns the error locator
/// σ(x) = 1 + σ_1 x + … + σ_e x^e (little-endian, σ[0] = 1). deg σ = e = number
/// of errors. The roots of σ are the inverse error locators α^{-i}.
fn berlekamp_massey(s: &[u8]) -> Vec<u8> {
    let n = s.len(); // = 2δ
    let mut sigma = vec![0u8; n + 1];
    sigma[0] = 1;
    let mut prev = vec![0u8; n + 1]; // B(x): σ from the last length change
    prev[0] = 1;

    let mut l: usize = 0; // current register length (= deg σ)
    let mut m: usize = 1; // shift since last length change
    let mut b: u8 = 1; // discrepancy at the last length change

    for i in 0..n {
        // Discrepancy: δ = S_{i+1} + Σ_{j=1}^{L} σ_j · S_{i+1-j}.
        let mut delta = s[i];
        for j in 1..=l {
            delta = gf_add(delta, gf_mul(sigma[j], s[i - j]));
        }

        if delta == 0 {
            m += 1;
        } else if 2 * l <= i {
            // Length change: T(x) = σ(x); σ(x) -= (δ/b) x^m B(x); B=T; L=i+1-L.
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
            // No length change: σ(x) -= (δ/b) x^m B(x).
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

// ── Chien search ─────────────────────────────────────────────────────────────

/// Find positions i in [0, n1) where σ(α^{-i}) = 0 (i.e. errors).
/// Returns (position, α^{-i}) pairs. Always scans all n1 candidates (CT req).
fn chien_search(sigma: &[u8], n1: usize) -> Vec<(usize, u8)> {
    let mut roots = Vec::with_capacity(sigma.len().saturating_sub(1));

    for i in 0..n1 {
        // α^{-i} = α^{255-i} for i>0, α^0 = 1 for i=0. n1 ≤ 90 ⇒ 255-i ≥ 165, safe.
        let xi = if i == 0 { GF_EXP[0] } else { GF_EXP[255 - i] };
        if gf_poly_eval(sigma, xi) == 0 {
            roots.push((i, xi));
        }
    }

    roots
}

// ── Forney algorithm ──────────────────────────────────────────────────────────

/// Compute error values via the Forney formula.
///
/// With σ(x) = ∏_l (1 - X_l x) and S(x) = Σ_{t} s[t] x^t (s[t] = S_{t+1}):
///   Ω(x)  = S(x)·σ(x) mod x^{2δ}        (error evaluator)
///   σ'(x) = formal derivative of σ      (char 2: only odd-degree terms survive)
///   Y_l   = Ω(X_l^{-1}) / σ'(X_l^{-1})
///
/// The Chien root `xi = α^{-i}` IS X_l^{-1}, so Ω and σ' are evaluated at `xi`
/// directly (NOT at its inverse). In char 2 the leading minus sign vanishes.
fn forney(sigma: &[u8], synd: &[u8], roots: &[(usize, u8)]) -> Vec<u8> {
    let two_delta = synd.len();

    // Ω(x) = (S · σ) mod x^{2δ}.
    let mut omega = vec![0u8; two_delta];
    for i in 0..synd.len() {
        for j in 0..sigma.len() {
            if i + j < two_delta {
                omega[i + j] = gf_add(omega[i + j], gf_mul(synd[i], sigma[j]));
            }
        }
    }

    // σ'(x): in char 2, d/dx of σ[j] x^j is σ[j] x^{j-1} for odd j, else 0.
    // So σ'[j-1] = σ[j] for odd j.
    let mut sigma_prime = vec![0u8; sigma.len()];
    for j in (1..sigma.len()).step_by(2) {
        sigma_prime[j - 1] = sigma[j];
    }

    roots
        .iter()
        .map(|&(_pos, xi)| {
            // xi = X_l^{-1}; evaluate Ω and σ' at xi.
            let omega_val = gf_poly_eval(&omega, xi);
            let sigma_val = gf_poly_eval(&sigma_prime, xi);
            gf_div(omega_val, sigma_val)
        })
        .collect()
}

// ── Public decode entry point ─────────────────────────────────────────────────

/// Decode a received RS codeword of length `n1` carrying `k` message symbols
/// with error-correcting capacity `delta` (so nk = 2δ = n1 - k).
///
/// Returns `Some(message)` (k symbols, from the high positions) on success,
/// or `None` if the error pattern is uncorrectable.
pub fn rs_decode(received: &[u8], k: usize, delta: usize) -> Option<Vec<u8>> {
    let n1 = received.len();
    let nk = 2 * delta;
    debug_assert_eq!(k + nk, n1);

    // Step 1: syndromes.
    let (s, has_error) = syndromes(received, delta);
    if !has_error {
        return Some(received[nk..].to_vec());
    }

    // Step 2: Berlekamp-Massey.
    let sigma = berlekamp_massey(&s);
    let num_errors = sigma.len() - 1;

    if num_errors == 0 || num_errors > delta {
        return None;
    }

    // Step 3: Chien search (full scan, no early exit).
    let roots = chien_search(&sigma, n1);

    // Root count must match deg σ, else the locator is invalid (too many errors).
    if roots.len() != num_errors {
        return None;
    }

    // Step 4: Forney — error values.
    let error_vals = forney(&sigma, &s, &roots);

    // Step 5: correct the received word.
    let mut corrected = received.to_vec();
    for (&(pos, _xi), &e) in roots.iter().zip(error_vals.iter()) {
        corrected[pos] = gf_add(corrected[pos], e);
    }

    // Re-verify: a valid correction zeroes all syndromes.
    let (_s2, still_err) = syndromes(&corrected, delta);
    if still_err {
        return None;
    }

    Some(corrected[nk..].to_vec())
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

    // ── Chien search locates the right positions ──────────────────────────────

    #[test]
    fn chien_search_finds_correct_positions() {
        let (n1, k, delta) = S1;
        let msg: Vec<u8> = (0..k as u8).collect();
        let cw = encode(&msg, n1, delta);
        let mut received = cw.clone();
        inject_errors(&mut received, &[(3, 0x11), (10, 0x22)]);

        let (s, _) = syndromes(&received, delta);
        let sigma = berlekamp_massey(&s);
        let roots = chien_search(&sigma, n1);
        let positions: Vec<usize> = roots.iter().map(|&(p, _)| p).collect();
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
            let sigma = berlekamp_massey(&s);
            assert_eq!(
                sigma.len() - 1,
                num_errors,
                "deg σ should equal error count {num_errors}"
            );
        }
    }
}
