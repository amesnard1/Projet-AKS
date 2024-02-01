#include<gmp.h>
#include<stdlib.h>
#include<math.h>
#include<stdio.h>

#include"small_primes.h"

unsigned long small_primes[23000] = tab_small_primes; // List of primes <= 2^18
unsigned long nb_small_primes = 23000;

// Integer power
unsigned long ipow(unsigned long r, unsigned long i) {
    if(i == 0)
        return 1;
    if(i % 2)
        return ipow(r, i - 1) * r;
    unsigned long p = ipow(r, i/2);
    return p * p;
}

// Computes log_2(n) +- epsilon.
// NOTE: assumes n > 2 and is not a power of 2
long double log2_mpz(mpz_t n, long double epsilon) {
    if(mpz_cmp_ui(n, 2) <= 0)
        return 1;
    unsigned long floor_log = mpz_sizeinbase(n, 2) - 1;
    long double lg = (long double) floor_log;
    mpf_t y_mpf, x_mpf;
    mpf_init(y_mpf);
    mpf_init(x_mpf);
    mpf_set_z(x_mpf, n);
    mpf_div_2exp(y_mpf, x_mpf, floor_log);
    mpf_t z_mpf;
    mpf_init(z_mpf);
    mpf_set(z_mpf, y_mpf);
    int m = 0;
    while(epsilon <= ((1.0) / (1 << (m+1)))) {
        while(mpf_cmp_ui(z_mpf, 2) < 0) {
            mpf_mul(z_mpf, z_mpf, z_mpf);
            m++;
        }
        if(m < 64)
            lg += pow((long double)2,(long double)-m);
        else
            return lg;
        mpf_div_2exp(z_mpf, z_mpf, 1);
    }
    return lg;
}

unsigned long log_square(mpz_t n) {
    long double lg = log2_mpz(n, 1e-10);
    lg += 1e-10; // Just to be sure we're slightly above log(n)
    lg *= lg;
    return (unsigned long) lg;
}

// Euler's totient
unsigned long euler(unsigned long r) {
    unsigned long e = 1;
    unsigned long i = 0;
    unsigned long p;
    unsigned long j;
    while ((r != 1) && (i < nb_small_primes)) {
        p = small_primes[i];
        if(r % p == 0) {
            j = 0;
            while(r % p == 0) {
                r /= p;
                j++;
            e *= p - 1;
            e *= ipow(p, j-1); } }
        i++; }
    if (r != 1)
      { e *= r - 1; }
    return e;
}

// This is the bound on a.
void bound(mpz_t target, unsigned long r, mpz_t n) {
    mpf_t bd;
    mpf_init(bd);
    mpf_sqrt_ui(bd, euler(r));
    mpf_t lg;
    mpf_init(lg);
    mpf_set_d(lg, log2_mpz(n, 1e-15));
    mpf_mul(bd, bd, lg);
    mpz_set_ui(target, mpf_get_ui(bd));
    mpf_clear(lg);
    mpf_clear(bd);
}

// Find the minimum r for which n has (multiplicative) order > log^2 n in Z_r.
// Assumes: n is not a perfect power. So either n = 1, 2, or log(n) is not an
// integer, so floor_log2(n) < log(n). Return 0 if n is found to be composite at
// this step, and 1 otherwise (an appropriate r is found)
unsigned long min_r(mpz_t n) {
    unsigned long lg = mpz_sizeinbase(n, 2) - 1;
    mpz_t d, lg2, k, power;
    mpz_init(d);
    mpz_init(lg2);
    mpz_init(k);
    mpz_init(power);
    for(unsigned long r = 2; r < (((unsigned long)1)<<36) ; r++) {
        if (mpz_cmp_ui(n ,r) <= 0) {
            mpz_clear(d);
            mpz_clear(lg2);
            mpz_clear(k);
            mpz_clear(power);
            return 1; }
        mpz_gcd_ui(d, n, r);
        if (mpz_cmp_ui(d, 1) !=0 ) {
            // mpz_cmp is zero if and only if the two members are equal
            mpz_clear(d);
            mpz_clear(lg2);
            mpz_clear(k);
            mpz_clear(power);
            return 0;
        }
        // Check if the order of n in Z_r is more than log^2.
        mpz_set_ui(lg2, mpz_sizeinbase(n, 2));
        mpz_mul(lg2, lg2, lg2);
        mpz_set_ui(power, 1);
        for(mpz_set_ui(k, 1); mpz_cmp_ui(k, log_square(n)) <= 0 && mpz_cmp_ui(k,
                    r) <= 0; mpz_add_ui(k, k, 1)) { 
            mpz_mul(power, power, n);
            mpz_mod_ui(power, power, r);
            if(mpz_cmp_ui(power, 1) == 0) {
                goto next_r; // This r doesn't work, skip to the next one (outer
                             // continue)
            }
        }
        // This r works, return it.
        mpz_clear(d);
        mpz_clear(lg2);
        mpz_clear(k);
        mpz_clear(power);
        return r;
        next_r:
    }
    // If r >= 2^36, we would need way too much memory
    printf("r is too big to handle\n");
    exit(1);
}
