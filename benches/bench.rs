// Criterion benchmarks for the performance-critical paths.
//
// Two groups:
//   • poly_mul — the ring multiplication hot path, both modes:
//       Mode A (mul_sparse_dense): sparse×dense, used in keygen / encrypt.
//       Mode B (mul_dense_ct):     constant-time dense×dense, used in decrypt.
//   • kem      — keygen / encaps / decaps end to end.
//
// Each benchmark runs for all three parameter sets (HQC-128 / 192 / 256), so
// the report shows how cost scales with the ring dimension n.
//
// Run with:  cargo bench
// A single group:  cargo bench --bench bench -- poly_mul
//
// These exist mainly to measure Step 17 (Karatsuba / SIMD poly_mul): capture a
// baseline now, then re-run after each optimization layer to confirm it pays
// off and to catch regressions.

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};

use hqcr::hash;
use hqcr::poly::mul::{mul_dense_ct, mul_sparse_dense};
use hqcr::poly::sampling::{sample_fixed_weight, sample_uniform};
use hqcr::poly::Poly;
use hqcr::kem::{decaps, encaps_deterministic, keygen_from_seed};
use hqcr::{Hqc128, Hqc192, Hqc256, HqcParams, SALT_BYTES, SEED_BYTES};

/// Build a deterministic dense (uniform) and sparse (weight-ω) operand pair for
/// the multiplication benchmarks, derived from a fixed seed so every run uses
/// identical inputs.
fn mul_operands<P: HqcParams>() -> (Poly<P>, Poly<P>) {
    let mut x = hash::xof(b"hqc-bench-poly-mul-seed");
    let dense = sample_uniform::<P>(&mut x);
    let sparse = sample_fixed_weight::<P>(&mut x, P::OMEGA);
    (dense, sparse)
}

fn bench_poly_mul(c: &mut Criterion) {
    let mut group = c.benchmark_group("poly_mul");

    fn run<P: HqcParams>(group: &mut criterion::BenchmarkGroup<criterion::measurement::WallTime>, name: &str) {
        let (dense, sparse) = mul_operands::<P>();

        // Mode A — sparse × dense (O(ω · N_WORDS)); keygen / encrypt path.
        group.bench_function(BenchmarkId::new("sparse_dense", name), |b| {
            b.iter(|| mul_sparse_dense::<P>(black_box(&sparse), black_box(&dense)))
        });

        // Mode B — constant-time dense × dense (O(N · N_WORDS)); decrypt path.
        group.bench_function(BenchmarkId::new("dense_ct", name), |b| {
            b.iter(|| mul_dense_ct::<P>(black_box(&dense), black_box(&sparse)))
        });
    }

    run::<Hqc128>(&mut group, "hqc128");
    run::<Hqc192>(&mut group, "hqc192");
    run::<Hqc256>(&mut group, "hqc256");

    group.finish();
}

fn bench_kem(c: &mut Criterion) {
    let mut group = c.benchmark_group("kem");

    fn run<P: HqcParams>(group: &mut criterion::BenchmarkGroup<criterion::measurement::WallTime>, name: &str) {
        let seed = [0x42u8; SEED_BYTES];
        let m = vec![0x11u8; P::K];
        let salt = [0x22u8; SALT_BYTES];

        // keygen_from_seed: SHAKE256 expansion + PKE.Keygen (s = x + h·y).
        group.bench_function(BenchmarkId::new("keygen", name), |b| {
            b.iter(|| keygen_from_seed::<P>(black_box(&seed)))
        });

        let (ek, dk) = keygen_from_seed::<P>(&seed);

        // encaps_deterministic: G + PKE.Encrypt (two Mode-A mults + C.Encode).
        group.bench_function(BenchmarkId::new("encaps", name), |b| {
            b.iter(|| encaps_deterministic::<P>(black_box(&ek), black_box(&m), black_box(&salt)))
        });

        // decaps: PKE.Decrypt (Mode-B mult + C.Decode) + re-encrypt + CT select.
        let (_k, ct) = encaps_deterministic::<P>(&ek, &m, &salt);
        group.bench_function(BenchmarkId::new("decaps", name), |b| {
            b.iter(|| decaps::<P>(black_box(&dk), black_box(&ct)))
        });
    }

    run::<Hqc128>(&mut group, "hqc128");
    run::<Hqc192>(&mut group, "hqc192");
    run::<Hqc256>(&mut group, "hqc256");

    group.finish();
}

criterion_group!(benches, bench_poly_mul, bench_kem);
criterion_main!(benches);
