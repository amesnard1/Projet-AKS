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

// Set to zero a polynomial of degree <= r-1
void zero_poly(unsigned long r, mpz_t* P) {
    for(unsigned long i = 0; i < r; i++)
        mpz_set_ui(P[i], 0);
}

// Copy polynomials of degree <= r-1
void copy_poly(unsigned long r, mpz_t* dest, const mpz_t* src) {
    for(unsigned long i = 0; i < r; i++)
        mpz_set(dest[i], src[i]);
}


// Multiply polynomials P and Q modulo X^r - 1 and store it in result.
// All the coefficients are seen modulo n.
// Assumes the polynomials are polynomials of length r
void naive_poly_mul_mod(unsigned long r, mpz_t n, mpz_t* result, unsigned long
        degP, const mpz_t* P, unsigned long degQ, const mpz_t* Q) {
    mpz_t tmp;
    mpz_init(tmp);
    for(unsigned long i = 0; i < degP + 1; i++) {
        for(unsigned long j = 0; j < degQ + 1; j++) {
            mpz_mul(tmp, P[i], Q[j]);
            mpz_mod(tmp, tmp, n);
            mpz_add(result[(i + j) % r], result[(i + j) % r], tmp);
            mpz_mod(result[(i + j) % r], result[(i + j) % r], n);
        }
    }
    mpz_clear(tmp);
}


// TODO use the naive multiplcation for small polynomials and implement FFT for
// larger polynomials. May start with Karatsuba, get a working implementation,
// then move to FFT to get the correct complexity, but perhaps it's a waste of
// time.

void poly_mul_mod(unsigned long r, mpz_t n, mpz_t* result, const mpz_t* P, const mpz_t* Q) {
    zero_poly(r, result);
    naive_poly_mul_mod(r, n, result, r-1, P, r-1, Q);
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
    for(long i = 0; i < r; i++)
        if(mpz_sgn(P[i]) != 0) return 0;
    return 1;
}

// Print the terms of degree <= k
// Yes, k is an unsigned long, but there's no human way to read anything that
// go beyond the bounds of an int.
void print_poly(unsigned long r, const mpz_t* P, unsigned long k) {
    int b = 0;
    if(P[0]) {
        gmp_printf("%Zd ", P[0]);
        b = 1;
    }
    if(P[1]) {
        if(b)
            printf("+ ");
        gmp_printf("%Zd x ", P[1]);
        b = 1;
    }
    for(unsigned long i = 2; i <= k; i++) {
        if(P[i]) {
            if(b)
                printf("+ ");
            gmp_printf("%Zd x^%d ", P[i], i);
            b = 1;
        }
    }
    printf("\n");
}

