/* main_speed.c -- KEM timing harness, shared by all six c_implementations folders.
 * Build from inside a folder: $(CC) $(CFLAGS) ../main_speed.c <objs> -I src -o bin/hqc-N-speed
 * Run: ./bin/hqc-N-speed [iterations]   (default 1000)
 *
 * Measures keygen / encaps / decaps wall-clock ns and rdtsc cycles.
 * Reports median, mean, min per operation.
 * Sizes every buffer from the CRYPTO_* macros in api.h (no hardcoded numbers).
 */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <x86intrin.h>  /* __rdtsc -- mingw-w64 gcc */

#include "api.h"
#include "symmetric.h"  /* prng_init */

/* ---- helpers ------------------------------------------------------------ */

static int cmp_uint64(const void *a, const void *b) {
    uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b;
    return (x > y) - (x < y);
}

static uint64_t median64(uint64_t *arr, int n) {
    qsort(arr, (size_t)n, sizeof *arr, cmp_uint64);
    return (n & 1) ? arr[n / 2] : (arr[n / 2 - 1] / 2 + arr[n / 2] / 2);
}

static uint64_t mean64(const uint64_t *arr, int n) {
    uint64_t s = 0;
    for (int i = 0; i < n; i++) s += arr[i];
    return s / (uint64_t)n;
}

/* wall-clock nanoseconds via CLOCK_MONOTONIC */
static uint64_t now_ns(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (uint64_t)ts.tv_sec * 1000000000ULL + (uint64_t)ts.tv_nsec;
}

static void print_stats(const char *label,
                        uint64_t *ns_arr, uint64_t *cy_arr, int n) {
    uint64_t med_ns = median64(ns_arr, n);
    uint64_t mn_ns  = ns_arr[0];           /* sorted after median64 */
    uint64_t avg_ns = mean64(ns_arr, n);
    uint64_t med_cy = median64(cy_arr, n);
    printf("  %-8s  med %8.3f ms  mean %8.3f ms  min %8.3f ms  med_cycles %10llu\n",
           label,
           (double)med_ns  / 1e6,
           (double)avg_ns  / 1e6,
           (double)mn_ns   / 1e6,
           (unsigned long long)med_cy);
}

/* ---- main --------------------------------------------------------------- */

int main(int argc, char *argv[]) {
    int N = (argc > 1) ? atoi(argv[1]) : 1000;
    if (N <= 0) N = 1000;

    /* deterministic seed -- mirrors the VERBOSE branch of main_hqc.c */
    uint8_t entropy[48];
    for (int i = 0; i < 48; i++) entropy[i] = (uint8_t)i;
    prng_init(entropy, NULL, sizeof entropy, 0);

    printf("=== %s  (N=%d iterations) ===\n", CRYPTO_ALGNAME, N);

    /* allocate buffers sized from api.h macros */
    uint8_t *pk  = malloc(CRYPTO_PUBLICKEYBYTES);
    uint8_t *sk  = malloc(CRYPTO_SECRETKEYBYTES);
    uint8_t *ct  = malloc(CRYPTO_CIPHERTEXTBYTES);
    uint8_t *ss1 = malloc(CRYPTO_BYTES);
    uint8_t *ss2 = malloc(CRYPTO_BYTES);
    if (!pk || !sk || !ct || !ss1 || !ss2) { fputs("malloc failed\n", stderr); return 1; }

    /* ---- correctness gate ---- */
    crypto_kem_keypair(pk, sk);
    crypto_kem_enc(ct, ss1, pk);
    crypto_kem_dec(ss2, ct, sk);
    if (memcmp(ss1, ss2, CRYPTO_BYTES) != 0) {
        printf("FAIL -- shared secrets do not match, aborting\n");
        return 1;
    }
    printf("PASS -- round-trip OK\n\n");

    /* anti-dead-code accumulator */
    volatile uint8_t acc = 0;

    uint64_t *ns = malloc((size_t)N * sizeof *ns);
    uint64_t *cy = malloc((size_t)N * sizeof *cy);
    if (!ns || !cy) { fputs("malloc failed\n", stderr); return 1; }

    /* ---- keygen ---- */
    for (int i = 0; i < N; i++) {
        uint64_t t0 = now_ns(), c0 = __rdtsc();
        crypto_kem_keypair(pk, sk);
        cy[i] = __rdtsc() - c0;
        ns[i] = now_ns() - t0;
        acc ^= pk[0];
    }
    print_stats("keygen", ns, cy, N);

    /* ---- encaps (pk fixed from last keygen) ---- */
    for (int i = 0; i < N; i++) {
        uint64_t t0 = now_ns(), c0 = __rdtsc();
        crypto_kem_enc(ct, ss1, pk);
        cy[i] = __rdtsc() - c0;
        ns[i] = now_ns() - t0;
        acc ^= ss1[0] ^ ct[0];
    }
    print_stats("encaps", ns, cy, N);

    /* ---- decaps (sk and ct fixed from above) ---- */
    for (int i = 0; i < N; i++) {
        uint64_t t0 = now_ns(), c0 = __rdtsc();
        crypto_kem_dec(ss2, ct, sk);
        cy[i] = __rdtsc() - c0;
        ns[i] = now_ns() - t0;
        acc ^= ss2[0];
    }
    print_stats("decaps", ns, cy, N);

    printf("\n(acc=%02x -- ignore, prevents dead-code elimination)\n", (unsigned)acc);

    free(pk); free(sk); free(ct); free(ss1); free(ss2); free(ns); free(cy);
    return 0;
}
