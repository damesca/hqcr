// HQC-KEM: IND-CCA2 KEM via salted FO with implicit rejection (SFO⊥_m).
//
// Keygen():
//   seed_KEM, σ ← random(32 B) each
//   (seed_dk, seed_ek) = SHA3-512(seed_KEM)
//   (ek, dk_PKE) = PKE.Keygen(seed_dk, seed_ek)
//   Return ek, dk = (ek, dk_PKE, σ, seed_KEM)
//
// Encaps(ek):
//   m ← random(k B), salt ← random(16 B)
//   (K, θ) = G(m, H(ek), salt)
//   c = PKE.Encrypt(ek, m, θ) || salt
//   Return (K, c)
//
// Decaps(dk, c):
//   m' = PKE.Decrypt(dk_PKE, c_PKE)        [may be None → use zero buffer for G]
//   (K', θ') = G(m' or zeros, H(ek), salt)
//   c' = PKE.Encrypt(ek, m', θ') || salt
//   K_bar = J(σ, c)
//   valid = CT: (!m'_is_none) & ct_eq(c', c)   [subtle::Choice + ConstantTimeEq]
//   Return ct_select(valid, K', K_bar)           [subtle::ConditionallySelectable]
