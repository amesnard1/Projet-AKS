#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

void poly_init(unsigned long r, mpz_t* P);
void poly_free(unsigned long r, mpz_t* P);
void zero_poly(unsigned long r, mpz_t* P);
void copy_poly(unsigned long r, mpz_t* dest, const mpz_t* src);
void naive_poly_mul_mod(unsigned long r, mpz_t n, mpz_t* result, unsigned long
        degP, const mpz_t* P, unsigned long degQ, const mpz_t* Q);
void karatsuba_poly_mul_mod(unsigned long r, mpz_t n, mpz_t* result, unsigned long
        degP, const mpz_t* P, unsigned long degQ, const mpz_t* Q, int square);
void poly_mul_mod(unsigned long r, mpz_t n, mpz_t* result, unsigned long n_P,
        const mpz_t* P, unsigned long n_Q, const mpz_t* Q, int square);
void poly_pow_mod(unsigned long r, mpz_t n, mpz_t* result, unsigned long n_P,
        const mpz_t* P, mpz_t k);
int is_zero_poly(unsigned long r, mpz_t* P);
void print_poly(unsigned long r, const mpz_t* P, unsigned long k);

#endif
