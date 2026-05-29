// Hash and XOF wrappers around the sha3 crate (RustCrypto).
// G  — SHAKE256: derives (θ, K) from (m, H(ek), salt); used in Encaps/Decaps.
// H  — SHA3-256: hashes ekKEM to a 32-byte digest used as input to G.
// J  — SHAKE256: derives the implicit rejection key from (σ, ciphertext).
// xof — SHAKE256 seed expander for sampling (returns impl XofReader).
// Note: SHA3-512 is used only for splitting seed_KEM into (seed_dk, seed_ek)
// in KEM keygen — it is NOT used for shared-key derivation.
