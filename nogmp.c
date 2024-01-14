#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include"arithmetic.h"
#include"polynomial.h"


// Check if the AKS equation for r, n, and a given.
int eqn(unsigned long r, mpz_t n, const mpz_t a) {
    mpz_t a_; mpz_init(a_);
    mpz_mod(a_, a, n);
    mpz_t n_mod_r_mpz;
    mpz_init(n_mod_r_mpz);
    mpz_mod_ui(n_mod_r_mpz, n, r);
    unsigned long n_mod_r = mpz_get_ui(n_mod_r_mpz);
    // The following two lines are essentially
    // why r has to be an unsigned long, not mpz_t
    mpz_t* P = (mpz_t*) malloc(r * sizeof(mpz_t));
    mpz_t* Q = (mpz_t*) malloc(r * sizeof(mpz_t));
    poly_init(r, P); poly_init(r, Q);
    zero_poly(r, P);
    zero_poly(r, Q);
    mpz_set(P[0], a);
    mpz_set_ui(P[1], 1);
    poly_pow_mod(r, n, Q, P, n);
    mpz_sub_ui(Q[n_mod_r], Q[n_mod_r], 1);
    mpz_sub(Q[0], Q[0], a);
    mpz_mod(Q[n_mod_r], Q[n_mod_r], n);
    mpz_mod(Q[0], Q[0], n);
    int b = is_zero_poly(r, Q);
    poly_free(r, P); poly_free(r, Q);
    mpz_clear(a_);
    return b;
}


// The actual AKS algorithm
int is_prime(mpz_t n) {
    if(mpz_perfect_power_p(n)) return 0; // AKS step 1
    unsigned long r = min_r(n);
    if(r == 0) return 0; // AKS steps 2 and 3
    if(mpz_cmp_ui(n, r) <= 0) return 1; // AKS step 4
    mpz_t limit, a, remainder;
    mpz_init(limit);
    mpz_init(a);
    mpz_init(remainder);
    bound(limit, r, n);
    for(mpz_set_ui(a, 1); mpz_cmp(a, limit) <= 0; mpz_add_ui(a, a, 1)) { // AKS step 5
        // TODO: add the possibility of a --verbose option to control logging
        //       levels.
        mpz_mod_ui(remainder, a, 50);
        if(!mpz_sgn(remainder))
            gmp_printf("progress: %Zd/%Zd\n", a, limit);
        if(!eqn(r, n, a)) return 0;
    }
    mpz_clear(limit);
    mpz_clear(a);
    mpz_clear(remainder);
    return 1;
}



int main() {
    mpz_t n;
    mpz_init(n);
    mpz_set_ui(n, 19997);
    mpz_t a;
    mpz_init(a);
    mpz_set_ui(a, 52);
    printf("%d\n", min_r(n));
    printf("%d\n", eqn(min_r(n), n, a));
    mpz_t i_mpz;
    mpz_init(i_mpz);
    int b;
    for(unsigned long int i = 1; i < 200; i++) {
        mpz_set_ui(i_mpz, i);
        b = is_prime(i_mpz);
        if(b)
            printf("n = %d; prime: %d\n", i, b);
    }
    mpz_set_ui(i_mpz, 672629);
    printf("%d", is_prime(i_mpz));
    mpz_clear(i_mpz);
    return 0;
}
