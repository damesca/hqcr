// RM(1,7) = [128, 8, 64] base code, duplicated P::MULTIPLICITY times (3× or 5×).
//
// Encoding: 8-bit input indexes the 128-bit row of the Hadamard matrix;
// codeword is then duplicated MULTIPLICITY times.
//
// Decoding (duplicated FHT):
//   1. Reshape into `multiplicity` sub-blocks of 128 bits.
//   2. Sum duplicates into i16 vector F of length 128 using (-1)^bit voting.
//   3. Apply length-128 Walsh-Hadamard Transform (WHT) over i16.
//   4. Decoded symbol = argmax |F̂|; sign determines the all-ones correction bit.
//   Tie-break: smallest index among positions with equal |F̂| (spec §3.4.3).
//
// CT requirement: argmax loop scans all 128 positions without early exit;
// running max updated via subtle::ConditionallySelectable (branchless).
