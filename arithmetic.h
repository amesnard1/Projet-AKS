#ifndef ARITHMETIC_H
#define ARITHMETIC_H

long floor_log2(long n);
long gcd(long a, long b);
long ipow(long a, long b);
long floor_root(long n, long b);
int is_perfect_power(long n);
long modpow(long r, long n, long k);
long modmul(long n, long a, long b);
long addmod(long n, long a, long b);
long bound(long r, long n);
long min_r(long n);

#endif
