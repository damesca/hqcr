// Shortened Reed-Solomon over GF(2^8).
// Generator polynomials G1, G2, G3 are hardcoded constants from spec §3.4.2
// (for RS-S1 [46,16,31], RS-S2 [56,24,33], RS-S3 [90,32,49] respectively).
//
// Encoding: systematic — message in first k positions, parity in remaining
// n1-k positions via polynomial long division over GF(2^8).
//
// Decoding:
//   1. Compute 2δ syndromes S_j = received(α^j) for j = 1..2δ.
//   2. Berlekamp-Massey → error locator polynomial σ(x).
//   3. Chien search: evaluate σ(α^{-i}) for all i in [0, n1) → error positions.
//   4. Forney algorithm → error values over GF(2^8).
//   Returns None if more than δ errors are detected.
//
// CT requirement: Chien search iterates all n1 candidates (46, 56, or 90)
// with no early exit, even after finding deg(σ) roots.
