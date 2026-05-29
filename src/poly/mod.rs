// Poly<P: HqcParams>: bit-packed polynomial in R = F2[X]/(X^n - 1).
// Stored as stack-allocated [u64; P::N_WORDS] (no heap).
// Operations: add (XOR), reduce (fold overflow bits), get_bit, set_bit, clear.
// Derives ZeroizeOnDrop for use as secret key material.
