#include<gmp.h>
#include<math.h>
#include<stdlib.h>
#include<stdio.h>

// TODO: Evaluate phi(r) and log n more rigourously for a tighter bound, if
// indeed found to be worthwhile to compute phi(r) (probably is, but needs
// testing).
void bound(mpz_t target, unsigned long r, mpz_t n) {
    mpz_set_ui(target, r - 1);
    mpz_sqrt(target, target);
    mpz_add_ui(target, target, 1); // approximately ceil(sqrt(r - 1))
    mpz_mul_ui(target, target, mpz_sizeinbase(n, 2));
}

// Find the minimum r for which n has (multiplicative) order > log^2 n in Z_r.
// Assumes: n is not a perfect power. So either n = 1, 2, or log(n) is not an
// integer, so floor_log2(n) < log(n). Return 0 if n is found to be composite at
// this step, and 1 otherwise (an appropriate r is found)
unsigned long min_r(mpz_t n) {
    unsigned long lg = mpz_sizeinbase(n, 2) - 1;
    mpz_t d, big_r, lg2, k, power;
    mpz_init(d);
    mpz_init(big_r);
    mpz_init(lg2);
    mpz_init(k);
    mpz_init(power);
    for(unsigned long r = 2; r > 0 ; r++) {
        mpz_set_ui(big_r, r);
        mpz_gcd(d, big_r, n);
        if(mpz_cmp_ui(d, 1) && mpz_cmp(d, n)) {
            // mpz_cmp is zero if and only if the two members are equal
            mpz_clear(d);
            mpz_clear(big_r);
            mpz_clear(lg2);
            mpz_clear(k);
            mpz_clear(power);
            return 0;
        }
        // Check if the order of n in Z_r is more than log^2.
        mpz_set_ui(lg2, mpz_sizeinbase(n, 2));
        mpz_mul(lg2, lg2, lg2);
        mpz_set(power, n);
        for(mpz_set_ui(k, 1); mpz_cmp(k, lg2) <= 0; mpz_add_ui(k, k, 1)) { 
            mpz_mul(power, power, n);
            mpz_mod(power, power, big_r);
            // mpz_powm(power, n, k, big_r);
            if(mpz_cmp_ui(power, 1) == 0) {
                goto next_r; // This r doesn't work, skip to the next one (outer
                             // continue)
            }
        }
        // This r works, return it.
        mpz_clear(d);
        mpz_clear(big_r);
        mpz_clear(lg2);
        mpz_clear(k);
        mpz_clear(power);
        return r;
        next_r:
    }
    // Here, r must have overflowed an unsigned int, which is a case we cannot
    // handle, so we will produce an error and exit the program.
    printf("r was not found before overflowing ulong, too big to handle\n");
    exit(1);
}
