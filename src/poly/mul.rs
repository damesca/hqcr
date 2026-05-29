// Polynomial multiplication in R — the hot path.
//
// Mode A (sparse × dense): used in keygen and encrypt; iterates over set-bit
// positions of the sparse operand and XOR-rotates the dense operand. O(ω·N/64).
//
// Mode B (dense × dense, CT): used in decrypt (u·y where y is secret); same
// as Mode A but with branchless bit extraction over all word positions to
// avoid leaking secret positions via timing.
//
// Optimization layers (in order):
//   L0 — portable word-level baseline (always compiled)
//   L1 — Karatsuba (~2× over L0 for Mode B)
//   L2 — SIMD pclmulqdq via std::arch::x86_64 (behind #[cfg(target_feature="pclmul")])
//
// unsafe is only permitted in the L2 SIMD block; L0 and L1 are safe Rust.
