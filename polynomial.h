#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

void zero_poly(long r, long* P);
void copy_poly(long r, long* dest, const long* src);
void naive_poly_mul_mod(long r, long n, long* result, long degP, const long* P, long degQ, const long* Q);
void poly_mul_mod(long r, long n, long* result, const long* P, const long* Q);
void poly_pow_mod(long r, long n, long* result, const long* P, long k);
int is_zero_poly(long r, long* P);
void print_poly(long r, const long* P, long k);

#endif
