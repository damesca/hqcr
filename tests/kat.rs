// KAT (Known Answer Test) integration tests.
// Parses official .rsp vector files from pqc-hqc.org and verifies byte-for-byte
// correctness of Keygen, Encaps, and Decaps for all three parameter sets
// (HQC-128, HQC-192, HQC-256). Enabled with `cargo test --features kat`.
// These are the sole correctness oracle — KAT passing means the implementation
// is correct per the spec.
