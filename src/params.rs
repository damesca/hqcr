// HqcParams sealed trait and three zero-size implementors: Hqc128, Hqc192, Hqc256.
// Defines all compile-time constants: N, N1, N2, K, OMEGA, OMEGA_R, DELTA,
// MULTIPLICITY, and derived sizes N_WORDS, PK_BYTES, CT_BYTES, SK_BYTES.

// Fixed sizes that do not vary across parameter sets.
pub const SEED_BYTES: usize = 32;
pub const SALT_BYTES: usize = 16;
pub const SHARED_KEY_BYTES: usize = 32;

// Sealed trait: only types in this crate can implement HqcParams.
mod sealed {
    pub trait Sealed {}
}

/// Compile-time parameters for a specific HQC security level.
///
/// The naming convention follows the 2025 spec:
///   N1 = external (RS) code length
///   N2 = internal (RM) code length after duplication
/// Some older documents swap these labels — always check against spec Table 5.
pub trait HqcParams: sealed::Sealed {
    /// Ring dimension n (a primitive prime so X^n - 1 has exactly two
    /// irreducible factors over F2, blocking algebraic attacks).
    const N: usize;

    /// Shortened RS code length (external code). Named n1 in the spec.
    const N1: usize;

    /// Duplicated RM code length (internal code). Named n2 in the spec.
    /// Equals 128 * MULTIPLICITY.
    const N2: usize;

    /// Message length in bytes.
    const K: usize;

    /// Weight of secret key vectors x and y.
    const OMEGA: usize;

    /// Weight of ephemeral vectors r1, r2, e (ωr = ωe in the spec).
    const OMEGA_R: usize;

    /// RS error-correcting capacity δ.
    const DELTA: usize;

    /// RM duplication factor: 3 for HQC-128, 5 for HQC-192 and HQC-256.
    const MULTIPLICITY: usize;

    /// Number of u64 words needed to store a polynomial of degree < N.
    /// ceil(N / 64).
    const N_WORDS: usize;

    /// Encapsulation key size in bytes: SEED_BYTES + ceil(N / 8).
    const PK_BYTES: usize;

    /// Ciphertext size in bytes: ceil(N/8) + ceil(N1*N2/8) + SALT_BYTES.
    const CT_BYTES: usize;

    /// Compressed secret key size (just seed_KEM): always 32 bytes.
    const SK_BYTES: usize;
}

// ── Parameter sets ────────────────────────────────────────────────────────────

/// HQC-128: NIST security level 1.
#[derive(Copy, Clone, Debug)]
pub struct Hqc128;

/// HQC-192: NIST security level 3.
#[derive(Copy, Clone, Debug)]
pub struct Hqc192;

/// HQC-256: NIST security level 5.
#[derive(Copy, Clone, Debug)]
pub struct Hqc256;

impl sealed::Sealed for Hqc128 {}
impl sealed::Sealed for Hqc192 {}
impl sealed::Sealed for Hqc256 {}

impl HqcParams for Hqc128 {
    const N: usize = 17_669;
    const N1: usize = 46;
    const N2: usize = 384; // 128 * 3
    const K: usize = 16;
    const OMEGA: usize = 66;
    const OMEGA_R: usize = 75;
    const DELTA: usize = 15;
    const MULTIPLICITY: usize = 3;

    const N_WORDS: usize = 17_669_usize.div_ceil(64);
    const PK_BYTES: usize = 32 + 17_669_usize.div_ceil(8);
    const CT_BYTES: usize = 17_669_usize.div_ceil(8) + ((46 * 384) as usize).div_ceil(8) + 16;
    const SK_BYTES: usize = 32;
}

impl HqcParams for Hqc192 {
    const N: usize = 35_851;
    const N1: usize = 56;
    const N2: usize = 640; // 128 * 5
    const K: usize = 24;
    const OMEGA: usize = 100;
    const OMEGA_R: usize = 114;
    const DELTA: usize = 16;
    const MULTIPLICITY: usize = 5;

    const N_WORDS: usize = 35_851_usize.div_ceil(64);
    const PK_BYTES: usize = 32 + 35_851_usize.div_ceil(8);
    const CT_BYTES: usize = 35_851_usize.div_ceil(8) + ((56 * 640) as usize).div_ceil(8) + 16;
    const SK_BYTES: usize = 32;
}

impl HqcParams for Hqc256 {
    const N: usize = 57_637;
    const N1: usize = 90;
    const N2: usize = 640; // 128 * 5
    const K: usize = 32;
    const OMEGA: usize = 131;
    const OMEGA_R: usize = 149;
    const DELTA: usize = 29;
    const MULTIPLICITY: usize = 5;

    const N_WORDS: usize = 57_637_usize.div_ceil(64);
    const PK_BYTES: usize = 32 + 57_637_usize.div_ceil(8);
    const CT_BYTES: usize = 57_637_usize.div_ceil(8) + ((90 * 640) as usize).div_ceil(8) + 16;
    const SK_BYTES: usize = 32;
}

// Reminder: ceil(a/b) = (a + b - a) / b

// ── Compile-time sanity checks ────────────────────────────────────────────────

const _: () = {
    assert!(Hqc128::PK_BYTES == 2_241);
    assert!(Hqc128::CT_BYTES == 4_433);
    assert!(Hqc128::N_WORDS == 277);

    assert!(Hqc192::PK_BYTES == 4_514);
    assert!(Hqc192::CT_BYTES == 8_978);
    assert!(Hqc192::N_WORDS == 561);

    assert!(Hqc256::PK_BYTES == 7_237);
    assert!(Hqc256::CT_BYTES == 14_421);
    assert!(Hqc256::N_WORDS == 901);
};
