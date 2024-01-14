#include<stdio.h>
#include<stdlib.h>
#include"arithmetic.h"

// Set to zero a polynomial of degree <= r-1
void zero_poly(long r, long* P) {
    for(long i = 0; i < r; i++)
        P[i] = 0;
}

// Copy polynomials of degree <= r-1
void copy_poly(long r, long* dest, const long* src) {
    for(long i = 0; i < r; i++)
        dest[i] = src[i];
}


// Multiply polynomials P and Q modulo X^r - 1 and store it in result.
// All the coefficients are seen modulo n.
// Assumes the polynomials are polynomials of length r
void naive_poly_mul_mod(long r, long n, long* result, long degP, const long* P, long degQ, const long* Q) {
    for(long i = 0; i < degP + 1; i++) {
        for(long j = 0; j < degQ + 1; j++) {
            result[(i + j) % r] = addmod(n, result[(i + j) % r], modmul(n, P[i], Q[j]));
        }
    }
}


// TODO use the naive multiplcation for small polynomials and implement FFT for
// larger polynomials. May start with Karatsuba, get a working implementation,
// then move to FFT to get the correct complexity, but perhaps it's a waste of
// time.

void poly_mul_mod(long r, long n, long* result, const long* P, const long* Q) {
    zero_poly(r, result);
    naive_poly_mul_mod(r, n, result, r-1, P, r-1, Q);
}

// Computes P^k mod X^r - 1 as a polynomial in Z_n[X]
void poly_pow_mod(long r, long n, long* result, const long* P, long k) {
    long* scratch = (long*) malloc(r * sizeof(long));
    long* base = (long*) malloc(r * sizeof(long));

    zero_poly(r, result);
    result[0] = 1;

    if(k == 0)
        return;

    for(long i = 0; i < r; i++) {
        scratch[i] = 0;
        base[i] = P[i];
    }

    while(k > 1) {
        if(k % 2 == 1) {
            poly_mul_mod(r, n, scratch, result, base);
            for(long i = 0; i < r; i++) {
                result[i] = scratch[i];
            }
            zero_poly(r, scratch);
            k--;
        }
        poly_mul_mod(r, n, scratch, base, base);
        for(long i = 0; i < r; i++) {
            base[i] = scratch[i];
        }
        zero_poly(r, scratch);
        k >>= 1;
    }

    poly_mul_mod(r, n, scratch, result, base);
    copy_poly(r, result, scratch);
    free(scratch);
    free(base);
}

// Check whether a polynomial of degree at most r-1 is the zero polynomial
int is_zero_poly(long r, long* P) {
    for(long i = 0; i < r; i++)
        if(P[i]) return 0;
    return 1;
}

// Print the terms of degree <= k
void print_poly(long r, const long* P, long k) {
    int b = 0;
    if(P[0]) {
        printf("%d ", P[0]);
        b = 1;
    }
    if(P[1]) {
        if(b)
            printf("+ ");
        printf("%d x ", P[1]);
        b = 1;
    }
    for(long i = 2; i <= k; i++) {
        if(P[i]) {
            if(b)
                printf("+ ");
            printf("%d x^%d ", P[i], i);
            b = 1;
        }
    }
    printf("\n");
}

