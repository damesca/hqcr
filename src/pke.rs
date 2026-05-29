// HQC-PKE: IND-CPA public-key encryption layer.
//
// Keygen(seed_dk, seed_ek):
//   Sample h = sample_uniform(SHAKE256(seed_ek))
//   Sample (y, x) = sample_fixed_weight × 2 from SHAKE256(seed_dk)  [y first]
//   Compute s = x + h·y
//   Return ek = (seed_ek, s), dk = (seed_dk, ek)
//
// Encrypt(ek, m, θ):
//   Sample (r2, e, r1) = sample_fixed_weight × 3 from SHAKE256(θ)  [r2 first]
//   Compute u = r1 + h·r2
//   Compute v = C.Encode(m) + s·r2 + e
//   Return (u, v)
//
// Decrypt(dk, (u, v)):
//   Compute tmp = v + u·y   [CT: Mode B multiplication]
//   Return C.Decode(tmp)    [None on decoding failure]
