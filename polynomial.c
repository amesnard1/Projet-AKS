#include<stdio.h>
#include<stdlib.h>
#include"arithmetic.h"

// Initialize all mpz_t entries of a polynomial
void poly_init(unsigned long r, mpz_t* P) {
    for(unsigned long i = 0; i < r; i++)
        mpz_init(P[i]);
}

// Clear all mpz_t entries of a polynomial and free the polynomial itself
void poly_free(unsigned long r, mpz_t* P) {
    for(unsigned long i = 0; i < r; i++)
        mpz_clear(P[i]);
    free(P);
}

// Set to zero a polynomial of degree <= nP - 1 (so nP is the number of elements in
// P)
void zero_poly(unsigned long nP,  mpz_t* P) {
    for(unsigned long i = 0; i < nP; i++)
        mpz_set_ui(P[i], 0);
}

// Copy polynomials of degree <= n_src - 1
void copy_poly(unsigned long n_src, mpz_t* dest, const mpz_t* src) {
    for(unsigned long i = 0; i < n_src; i++)
        mpz_set(dest[i], src[i]);
}

// Computes, in Z_n[X], P + Q X^shift modulo X^r - 1. Assumes result has enough
// allocated space // (at least max(n_P, n_Q + shift)).
void shifted_add_mod_poly(unsigned long r, mpz_t n, mpz_t* result, unsigned long
        n_P, const mpz_t* P, unsigned long n_Q, const mpz_t* Q, unsigned long
        shift) {
    mpz_t mpz_zero;
    mpz_init(mpz_zero);
    mpz_set_ui(mpz_zero, 0); // I think mpz_init does that anyway
    for(unsigned long i = 0; i < (n_P < n_Q ? n_P : n_Q); i++) {
        mpz_add(result[(i + shift)%r], (i < n_P) ? P[i] : mpz_zero, (i < n_Q) ?
                Q[i] : mpz_zero);
        mpz_mod(result[(i + shift)%r], result[i], n);
    }
    mpz_clear(mpz_zero);
}

// Add polynomials in Z_n[X]. Assumes result has enough allocated space (at
// least max(n_P, n_Q)).
void add_mod_poly(mpz_t n, mpz_t* result, unsigned long n_P, const mpz_t* P,
        unsigned long n_Q, const mpz_t* Q) {
    shifted_add_mod_poly(((n_Q < n_P) ? n_P : n_Q) + 5, n, result, n_P, P, n_Q, Q, 0);
}


// Subtract polynomials in Z_n[X] (does P - Q). Assumes result has enough
// allocated space (at least max(n_P, n_Q)).
void sub_mod_poly(mpz_t n, mpz_t* result, unsigned long n_P, const mpz_t* P,
        unsigned long n_Q, const mpz_t* Q) {
    mpz_t mpz_zero;
    mpz_init(mpz_zero);
    mpz_set_ui(mpz_zero, 0); // I think mpz_init does that anyway
    for(unsigned long i = 0; i < (n_P < n_Q ? n_P : n_Q); i++) {
        mpz_sub(result[i], (i < n_P) ? P[i] : mpz_zero, (i < n_Q) ? Q[i] : mpz_zero);
        mpz_mod(result[i], result[i], n);
    }
    mpz_clear(mpz_zero);
}


// Multiply polynomials P and Q modulo X^r - 1 and store it in result.
// All the coefficients are seen modulo n.
// Assumes the polynomials are polynomials of length r
// n_P = deg P + 1
void naive_poly_mul_mod(unsigned long r, mpz_t n, mpz_t* result, unsigned long
        n_P, const mpz_t* P, unsigned long n_Q, const mpz_t* Q) {
    mpz_t tmp;
    mpz_init(tmp);
    for(unsigned long i = 0; i < n_P; i++) {
        for(unsigned long j = 0; j < n_Q; j++) {
            mpz_mul(tmp, P[i], Q[j]);
            mpz_mod(tmp, tmp, n);
            mpz_add(result[(i + j) % r], result[(i + j) % r], tmp);
            mpz_mod(result[(i + j) % r], result[(i + j) % r], n);
        }
    }
    mpz_clear(tmp);
}


void karatsuba_poly_mul_mod(unsigned long r, mpz_t n, mpz_t* result, unsigned
        long n_P, const mpz_t* P, unsigned long n_Q, const mpz_t* Q) {
    // If the polynomials are small enough, favor the trivial multiplication
    // which would be faster.
    if(n_P < 30 || n_Q < 30) {
        naive_poly_mul_mod(r, n, result, n_P, P, n_Q, Q);
        return;
    }
    // All the polynomials involved here have degree at most degP + degQ
    unsigned long max_n_poly = n_P + n_Q;
    if(max_n_poly > r) max_n_poly = r;
    unsigned long m = (n_P < n_Q ? n_P : n_Q)/ 2;
    // P = P0 + X^m P1
    // Q = Q0 + X^m Q1
    mpz_t* P0 = P; // Same pointer, but will be treated as an array of m entries
                   // instead of r.
    mpz_t* P1 = P + m; // Essentially divides by X^m by shifting.
    mpz_t* Q0 = Q;
    mpz_t* Q1 = Q + m;

    mpz_t* P0Q0 = (mpz_t*) malloc(max_n_poly * sizeof(mpz_t));
    mpz_t* P1Q1 = (mpz_t*) malloc(max_n_poly * sizeof(mpz_t));
    // R = (P1 + P0) (Q1 + Q0)
    mpz_t* R = (mpz_t*) malloc(max_n_poly * sizeof(mpz_t));
    // RP = (P1 + P0)
    mpz_t* RP = (mpz_t*) malloc(max_n_poly * sizeof(mpz_t));
    // RQ = (Q1 + Q0)
    mpz_t* RQ = (mpz_t*) malloc(max_n_poly * sizeof(mpz_t));

    poly_init(max_n_poly, P0Q0);
    poly_init(max_n_poly, P1Q1);
    poly_init(max_n_poly, R);
    poly_init(max_n_poly, RP);
    poly_init(max_n_poly, RQ);

    zero_poly(max_n_poly, P0Q0);
    zero_poly(max_n_poly, P1Q1);
    zero_poly(max_n_poly, R);
    zero_poly(max_n_poly, RP);
    zero_poly(max_n_poly, RQ);

    karatsuba_poly_mul_mod(r, n, P0Q0, m, P0, m, Q0);
    karatsuba_poly_mul_mod(r, n, P1Q1, n_P - m, P1, n_P - m, Q1);
    add_mod_poly(n, RP, m, P0, n_P - m, P1);
    add_mod_poly(n, RQ, m, Q0, n_Q - m, Q1);
    karatsuba_poly_mul_mod(r, n, R, (m < n_P - m ? n_P - m : m), RP,
                                    (m < n_Q - m ? n_Q - m : m), RQ);
    // R = R - P0Q0 - P1Q1
    sub_mod_poly(n, R, max_n_poly, R, 2 * m, P0Q0);
    sub_mod_poly(n, R, max_n_poly, R, (n_P - m) + (n_Q - m), P0Q0);

    zero_poly(max_n_poly, result);
    shifted_add_mod_poly(r, n, result, max_n_poly, result, max_n_poly, P0Q0, 0);
    shifted_add_mod_poly(r, n, result, max_n_poly, result, max_n_poly, R, m);
    shifted_add_mod_poly(r, n, result, max_n_poly, result, max_n_poly, P1Q1, 2*m);
    poly_free(max_n_poly, P0Q0);
    poly_free(max_n_poly, P1Q1);
    poly_free(max_n_poly, R);
    poly_free(max_n_poly, RP);
    poly_free(max_n_poly, RQ);
}


// TODO use the naive multiplcation for small polynomials and implement FFT for
// larger polynomials. May start with Karatsuba, get a working implementation,
// then move to FFT to get the correct complexity, but perhaps it's a waste of
// time.
void poly_mul_mod(unsigned long r, mpz_t n, mpz_t* result, const mpz_t* P, const mpz_t* Q) {
    zero_poly(r, result);
    karatsuba_poly_mul_mod(r, n, result, r, P, r, Q);
}

// Computes P^k mod X^r - 1 as a polynomial in Z_n[X] using an iterative
// square and multiply algorithm.
void poly_pow_mod(unsigned long r, mpz_t n, mpz_t* result, const mpz_t* P, mpz_t k) {
    mpz_t* scratch = (mpz_t*) malloc(r * sizeof(mpz_t));
    mpz_t* base = (mpz_t*) malloc(r * sizeof(mpz_t));

    poly_init(r, scratch);
    poly_init(r, base);

    zero_poly(r, result);
    mpz_set_ui(result[0], 1);

    if(mpz_cmp_ui(k, 0) == 0)
        return;

    for(unsigned long i = 0; i < r; i++) {
        mpz_set_ui(scratch[i], 0);
        mpz_set(base[i], P[i]);
    }

    while(mpz_cmp_ui(k, 1) > 0) {
        if(mpz_odd_p(k)) {
            poly_mul_mod(r, n, scratch, result, base);
            copy_poly(r, result, scratch);
            zero_poly(r, scratch);
            mpz_sub_ui(k, k, 1);
        }
        poly_mul_mod(r, n, scratch, base, base);
        copy_poly(r, base, scratch);
        zero_poly(r, scratch);
        mpz_fdiv_q_2exp(k, k, 1);
    }

    poly_mul_mod(r, n, scratch, result, base);
    copy_poly(r, result, scratch);
    poly_free(r, scratch);
    poly_free(r, base);
}

// Check whether a polynomial of degree at most r-1 is the zero polynomial
int is_zero_poly(unsigned long r, mpz_t* P) {
    for(unsigned long i = 0; i < r; i++)
        if(mpz_sgn(P[i]) != 0) return 0;
    return 1;
}

// Print the terms of degree <= k
// Yes, k is an unsigned long, but there's no human way to read anything that
// go beyond the bounds of an int.
//void print_poly(unsigned long r, const mpz_t* P, unsigned long k) {
//    int b = 0;
//    if(P[0]) {
//        gmp_printf("%Zd ", P[0]);
//        b = 1;
//    }
//    if(P[1]) {
//        if(b)
//            printf("+ ");
//        gmp_printf("%Zd x ", P[1]);
//        b = 1;
//    }
//    for(unsigned long i = 2; i <= k; i++) {
//        if(P[i]) {
//            if(b)
//                printf("+ ");
//            gmp_printf("%Zd x^%d ", P[i], i);
//            b = 1;
//        }
//    }
//    printf("\n");
//}

