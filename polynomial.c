#include<stdio.h>
#include<stdlib.h>
#include"arithmetic.h"

// Print polynomial
void print_poly(unsigned long r, const mpz_t* P) {
    int print_started = 0;
    if(r <= 0) goto end;
    if(mpz_sgn(P[0]) != 0) {
        print_started = 1;
        gmp_printf("%Zd", P[0]);
    }
    if(r <= 1) goto end;
    if(mpz_sgn(P[1]) != 0) {
        if(print_started)
            printf(" + ");
        print_started = 1;
        gmp_printf("%Zd X", P[1]);
    }
    for(unsigned long i = 2; i < r; i++) {
        if(mpz_sgn(P[i]) == 0) continue;
        if(print_started)
            printf(" + ");
        print_started = 1;
        gmp_printf("%Zd X^%d", P[i], i);
    }
end:
    if(!print_started)
        printf("0");
    printf("\n");
}
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
// allocated space (at least max(n_P, n_Q + shift)).
//
// !!!!!!!!!!!
// NOTE: This function does not zero the vector result and will add on top of
// it, meaning that it will effectively do the operation
// result = result + P + Q X^shift.
//
// In particular, if it is used as add(result, result, P) with the intention of
// meaning: result = result + P, it will actually do the operation:
// result = result + result + P = 2 * result + P.
//
// For this exact usecase, a boolean ignore_P can be used to instead evaluate
// result = result + Q. As a safeguard, ignore_P will be ignored unless result =
// P (as pointers)
//
// Tracking the length of result for zeroing appears to be too annoying in this
// form. Since size tracking is so common across this project, it may be
// worthwhile to introduce a polynomial struct, or perhaps even a class.
// !!!!!!!!!!!
void shifted_add_mod_poly(unsigned long r, const mpz_t n, mpz_t* result, unsigned long
        n_P, const mpz_t* P, unsigned long n_Q, const mpz_t* Q, unsigned long
        shift, int ignore_P) {
    mpz_t mpz_zero;
    mpz_init(mpz_zero);
    mpz_set_ui(mpz_zero, 0); // I think mpz_init does that anyway
    for(unsigned long i = 0; i < (n_P < n_Q + shift ? n_Q + shift: n_P); i++) {
        if(!ignore_P || P != result)
            mpz_add(result[i%r], result[i%r], (i < n_P) ? P[i] : mpz_zero);
        mpz_add(result[i%r], result[i%r], (shift <= i && i < n_Q + shift) ?  Q[i
                - shift] : mpz_zero);
        mpz_mod(result[i%r], result[i%r], n);
    }
    mpz_clear(mpz_zero);
}

// Add polynomials in Z_n[X]. Assumes result has enough allocated space (at
// least max(n_P, n_Q)). Check shifted_add_mod_poly for a definition of
// ignore_P.
void add_mod_poly(unsigned long r, const mpz_t n, mpz_t* result, unsigned long
        n_P, const mpz_t* P, unsigned long n_Q, const mpz_t* Q, int ignore_P) {
    shifted_add_mod_poly(r, n, result, n_P, P, n_Q, Q, 0, ignore_P);
}


// Subtract polynomials in Z_n[X] (does P - Q). Assumes result has enough
// allocated space (at least max(n_P, n_Q)).
void sub_mod_poly(const mpz_t n, mpz_t* result, unsigned long n_P, const mpz_t* P,
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


// Evaluates Karatsuba multiplication on polynomials, with a special case for
// squaring (same asymptotic complexity, but less work copying)
void karatsuba_poly_mul_mod(unsigned long r, const mpz_t n, mpz_t* result, unsigned
        long n_P, const mpz_t* P, unsigned long n_Q, const mpz_t* Q, int square) {
    // If the polynomials are small enough, favor the trivial multiplication
    // which would be faster.
    if(n_P < 40 || n_Q < 40) {
        naive_poly_mul_mod(r, n, result, n_P, P, n_Q, Q);
        return;
    }
    // All the polynomials involved here have degree at most degP + degQ
    unsigned long max_n_poly = n_P + n_Q;
    if(max_n_poly > r) max_n_poly = r;
    unsigned long m = (n_P < n_Q ? n_P : n_Q)/ 2;
    if(m > r/2) m = r/2;
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
    if(!square) {
        poly_init(max_n_poly, RP);
        poly_init(max_n_poly, RQ);
    }

    karatsuba_poly_mul_mod(r, n, P0Q0, m, P0, m, Q0, square);
    karatsuba_poly_mul_mod(r, n, P1Q1, n_P - m, P1, n_Q - m, Q1, square);
    if(square) {
        karatsuba_poly_mul_mod(r, n, R, m, P0, n_P - m, P1, 0);
        add_mod_poly(r, n, R, max_n_poly, R, max_n_poly, R, 1);
    }
    else {
        add_mod_poly(r, n, RP, m, P0, n_P - m, P1, 0);
        add_mod_poly(r, n, RQ, m, Q0, n_Q - m, Q1, 0);
        karatsuba_poly_mul_mod(r, n, R, (m < n_P - m ? n_P - m : m), RP,
                                        (m < n_Q - m ? n_Q - m : m), RQ, 0);
        // R = R - P0Q0 - P1Q1
        sub_mod_poly(n, R, max_n_poly, R, 2 * m, P0Q0);
        sub_mod_poly(n, R, max_n_poly, R, (n_P - m) + (n_Q - m), P1Q1);
    }

    zero_poly(max_n_poly, result);
    shifted_add_mod_poly(r, n, result, max_n_poly, result, max_n_poly, P0Q0, 0, 1);
    shifted_add_mod_poly(r, n, result, max_n_poly, result, max_n_poly, R, m, 1);
    shifted_add_mod_poly(r, n, result, max_n_poly, result, max_n_poly, P1Q1, 2*m, 1);

    poly_free(max_n_poly, P0Q0);
    poly_free(max_n_poly, P1Q1);
    poly_free(max_n_poly, R);
    if(!square) {
        poly_free(max_n_poly, RP);
        poly_free(max_n_poly, RQ);
    }
}


// Multiply polynomials P and Q and store into result.
void poly_mul_mod(unsigned long r, mpz_t n, mpz_t* result, unsigned long n_P,
        const mpz_t* P, unsigned long n_Q, const mpz_t* Q, int square) {
    zero_poly(r, result);
    karatsuba_poly_mul_mod(r, n, result, n_P, P, n_Q, Q, square);
}

// Computes P^k mod X^r - 1 as a polynomial in Z_n[X] using an iterative
// square and multiply algorithm.
void poly_pow_mod(unsigned long r, mpz_t n, mpz_t* result, unsigned long n_P,
        const mpz_t* P, const mpz_t k) {
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

    unsigned long n_base = n_P;
    unsigned long n_result = 1;

    mpz_t k_;
    mpz_init(k_);
    mpz_set(k_, k);
    while(mpz_cmp_ui(k_, 1) > 0) {
        if(mpz_odd_p(k_)) {
            poly_mul_mod(r, n, scratch, n_result, result, n_base, base, 0);
            copy_poly(r, result, scratch);
            n_result = (n_base + n_result);
            n_result = (n_result > r) ? r : n_result;
            zero_poly(r, scratch);
            mpz_sub_ui(k_, k_, 1);
        }
        poly_mul_mod(r, n, scratch, n_base, base, n_base, base, 1);
        copy_poly(r, base, scratch);
        n_base = 2 * n_base;
        n_base = (n_base > r) ? r : n_base;
        zero_poly(r, scratch);
        mpz_fdiv_q_2exp(k_, k_, 1);
    }

    mpz_clear(k_);
    poly_mul_mod(r, n, scratch, n_result, result, n_base, base, 0);
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

