#ifndef ARITHMETIC_H
#define ARITHMETIC_H

#include<gmp.h>

// unsigned long is used for r because we will need to store r-long arrays in
// memory, so it would be pointless to go beyond unsigned long since there are
// no reasonable memories that large.

void bound(mpz_t target, unsigned long r, mpz_t n);
unsigned long min_r(mpz_t n);

#endif
