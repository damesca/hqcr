//! # hqc — Hamming Quasi-Cyclic post-quantum KEM
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
//! [`encaps`], and [`decaps`] (plus the deterministic variants
//! [`keygen_from_seed`] and [`encaps_deterministic`]). The shared secret is a
//! [`SharedKey`] (`[u8; 32]`).
//!
//! With your own CSPRNG (anything implementing `rand_core::{RngCore, CryptoRng}`,
//! e.g. `rand::rngs::OsRng`):
//!
//! ```no_run
//! use hqc::Hqc128;
//!
//! # fn demo<R: rand_core::RngCore + rand_core::CryptoRng>(rng: &mut R) {
//! let (ek, dk) = hqc::keygen::<Hqc128, _>(rng);     // public + secret key
//! let (k_send, ct) = hqc::encaps::<Hqc128, _>(rng, &ek); // shared secret + capsule
//! let k_recv = hqc::decaps::<Hqc128>(&dk, &ct);     // recovered shared secret
//! assert_eq!(k_send, k_recv);
//! # }
//! ```
//!
//! The deterministic API takes the randomness explicitly — reproducible, and
//! used by the KAT harness. This example actually runs and round-trips:
//!
//! ```
//! use hqc::Hqc128;
//!
//! // In production, draw `seed`, `m`, and `salt` from a CSPRNG (see `keygen`).
//! let seed = [0x42u8; hqc::SEED_BYTES];
//! let (ek, dk) = hqc::keygen_from_seed::<Hqc128>(&seed);
//!
//! let m = [0x11u8; 16]; // Hqc128::K == 16 bytes
//! let salt = [0x22u8; hqc::SALT_BYTES];
//! let (k_send, ct) = hqc::encaps_deterministic::<Hqc128>(&ek, &m, &salt);
//!
//! let k_recv = hqc::decaps::<Hqc128>(&dk, &ct);
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
// `hqc::keygen` / `hqc::encaps` / `hqc::decaps` rather than reaching into
// `hqc::kem::`.
pub use kem::{DecapsulationKey, PublicKey};
pub use kem::{decaps, encaps, encaps_deterministic, keygen, keygen_from_seed};
