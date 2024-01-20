#ifndef ARITHMETIC_H
#define ARITHMETIC_H

#include<gmp.h>

// unsigned long is used for r because we will need to store r-long arrays in
// memory, so it would be pointless to go beyond unsigned long since there are
// no reasonable memories that large.

unsigned long ipow(unsigned long, unsigned long);
long double log2_mpz(mpz_t n, long double epsilon);
unsigned long euler(unsigned long);
void bound(mpz_t target, unsigned long r, mpz_t n);
unsigned long min_r(mpz_t n);

#endif
