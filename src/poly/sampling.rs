// Sampling functions for polynomials in R.
//
// sample_fixed_weight: rejection sampling producing exactly `weight` distinct
// positions in [0, n). CT requirement: deduplication check uses
// subtle::ConstantTimeEq; branching only on the public condition pos < n.
//
// sample_uniform: fills all N bits from a SHAKE256 XOF — used for public h.
//
// hamming_weight: popcount over u64 limbs via u64::count_ones().
