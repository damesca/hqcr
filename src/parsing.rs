// Byte-level serialization matching the spec wire format exactly (validated by KAT).
// Bit-pack/unpack Poly<P> to/from [u8; ceil(N/8)].
// Public key layout:  seed_h (32 B) || s (ceil(n/8) B)
// Ciphertext layout:  u (ceil(n/8) B) || v (ceil(n1·128/8) B) || salt (16 B)
// Secret key (compressed): seed_KEM (32 B)
// Secret key (full):  ekKEM || dkPKE || σ (32 B) || seed_KEM (32 B)
// Trailing bits of v (ℓ = n - n1·n2 bits) are always zeroed on serialization.
