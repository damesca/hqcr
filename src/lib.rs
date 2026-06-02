// Crate root: re-exports of public API, HqcParams trait definition (sealed),
// and top-level KEM entry points for Hqc128, Hqc192, Hqc256.

pub mod params;
pub(crate) mod gf;
pub mod poly;

pub use params::{HqcParams, Hqc128, Hqc192, Hqc256};
pub use params::{SEED_BYTES, SALT_BYTES, SHARED_KEY_BYTES};
