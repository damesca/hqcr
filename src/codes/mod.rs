// Top-level concatenated codec: C.Encode and C.Decode.
// Encoding order: RS (outer) then RM (inner) — RS.Encode followed by RM.Encode
// for each symbol, producing n1 * (n2/n1) bits embedded in a Poly<P>.
// Decoding order: RM (inner) then RS (outer) — split into n1 blocks, RM.Decode
// each, then RS.Decode the resulting n1 GF(2^8) symbols back to k bytes.
