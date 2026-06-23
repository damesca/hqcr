//! # hqcr — Hamming Quasi-Cyclic post-quantum KEM
//!
//! A pure-Rust implementation of **HQC**, the code-based key-encapsulation
//! mechanism selected by NIST in March 2025. It provides the IND-CCA2 KEM
//! (built on an IND-CPA PKE via the salted Fujisaki–Okamoto transform with
//! implicit rejection) for all three parameter sets — [`Hqc128`], [`Hqc192`],
//! and [`Hqc256`].
//!
//! ## ⚠️ Not production-ready
//!
//! This crate is a **learning project**. It is validated byte-for-byte against
//! the official NIST KAT vectors, but it has **not** had an independent security
//! review or side-channel audit, and was developed with AI assistance. Do not
//! use it to protect anything real.
//!
//! ## Quick start
//!
//! The high-level entry points are re-exported at the crate root: [`keygen`],
//! [`encaps`], and [`decaps`]. The shared secret is a
//! [`SharedKey`] (`[u8; 32]`).
//!
//! With your own CSPRNG (anything implementing `rand_core::{RngCore, CryptoRng}`,
//! e.g. `rand::rngs::OsRng`):
//!
//! ```no_run
//! use hqcr::Hqc128;
//!
//! # fn demo<R: rand_core::RngCore + rand_core::CryptoRng>(rng: &mut R) {
//! let (ek, dk) = hqcr::keygen::<Hqc128, _>(rng);     // public + secret key
//! let (k_send, ct) = hqcr::encaps::<Hqc128, _>(rng, &ek); // shared secret + capsule
//! let k_recv = hqcr::decaps::<Hqc128>(&dk, &ct);     // recovered shared secret
//! assert_eq!(k_send, k_recv);
//! # }
//! ```
//!
//! The deterministic API takes the randomness explicitly — reproducible, and
//! used by the KAT harness. This example actually runs and round-trips:
//!
//! ```
//! use hqcr::Hqc128;
//!
//! // In production, draw `seed`, `m`, and `salt` from a CSPRNG (see `keygen`).
//! let seed = [0x42u8; hqcr::SEED_BYTES];
//! let (ek, dk) = hqcr::kem::keygen_from_seed::<Hqc128>(&seed);
//!
//! let m = [0x11u8; 16]; // Hqc128::K == 16 bytes
//! let salt = [0x22u8; hqcr::SALT_BYTES];
//! let (k_send, ct) = hqcr::kem::encaps_deterministic::<Hqc128>(&ek, &m, &salt);
//!
//! let k_recv = hqcr::decaps::<Hqc128>(&dk, &ct);
//! assert_eq!(k_send, k_recv);
//! ```
//!
//! ## Module map
//!
//! [`kem`] is the public surface for almost all callers. [`pke`] exposes the
//! lower-level IND-CPA layer (**not** CCA-secure on its own). [`params`] holds
//! the [`HqcParams`] trait and the three parameter markers. The remaining
//! modules ([`poly`], [`codes`], [`hash`], [`parsing`]) are implementation
//! internals — `pub` so the test harnesses can reach them, but hidden from these
//! docs and not part of the stable API.

pub mod params;
pub(crate) mod gf;
pub mod pke;
pub mod kem;

// Implementation internals: kept `pub` so the integration / KAT / intermediate
// harnesses in `tests/` can reach them, but hidden from the rendered docs and
// not covered by any API-stability promise.
#[doc(hidden)]
pub mod poly;
#[doc(hidden)]
pub mod codes;
#[doc(hidden)]
pub mod parsing;
#[doc(hidden)]
pub mod hash;

// ── Public API surface ────────────────────────────────────────────────────────

pub use params::{HqcParams, Hqc128, Hqc192, Hqc256};
pub use params::{SALT_BYTES, SEED_BYTES, SHARED_KEY_BYTES};

/// The 32-byte shared secret produced by [`encaps`] and recovered by [`decaps`].
pub use hash::SharedKey;

// KEM key types and entry points, flattened to the crate root so callers write
// `hqcr::keygen` / `hqcr::encaps` / `hqcr::decaps` rather than reaching into
// `hqcr::kem::`.
pub use kem::{DecapsulationKey, PublicKey};
pub use kem::{decaps, encaps, keygen};

// ── Zeroize audit: compile-time guard (19b, Layer 2) ──────────────────────────
//
// Trip-wire for the zeroize invariant (docs/audit/zeroize.md §5): if a secret
// type ever loses its `#[derive(ZeroizeOnDrop)]`, the bound below fails to
// compile. And because the derived `ZeroizeOnDrop` brings a `Drop` impl — which
// Rust forbids alongside `Copy` — this also guarantees these types can never
// silently become `Copy` (the bitwise-copy defeater). Test-only; no effect on a
// normal build.
#[cfg(test)]
mod zeroize_guards {
    use crate::params::{Hqc128, Hqc192, Hqc256};
    use crate::pke::DecryptionKey;
    use crate::poly::Poly;

    fn assert_zeroize_on_drop<T: zeroize::ZeroizeOnDrop>() {}

    #[test]
    fn secret_types_are_zeroize_on_drop() {
        assert_zeroize_on_drop::<Poly<Hqc128>>();
        assert_zeroize_on_drop::<Poly<Hqc192>>();
        assert_zeroize_on_drop::<Poly<Hqc256>>();
        assert_zeroize_on_drop::<DecryptionKey>();
    }
}
