#include<math.h>

// Returns the floor of the log base 2 of n if n >= 1, -1 if n == 0
// Complexity: log_2 n
long floor_log2(long n) {
  long lg = -1;
  while(n) {
      lg++;
      n >>= 1;
  }
  return lg;
}

// Replace with GMP GCD
long gcd(long a, long b) {
    while (b != 0) {
        long temp = b;
        b = a % b;
        a = temp;
    }
    return a;
}

// Inspired by Wikipedia ("Exponentiation by squaring") because for some reason
// I kept on getting lost. This iterative structure for square and multiply is
// important for polynomial powers later.
long ipow(long a, long b) {
    if(b == 0) return 1;
    long result = 1;
    while(b > 1) {
        if(b % 2 == 1) {
            result = result * a;
            b--;
        }
        a = a * a;
        b >>= 1;
    }
    return result * a;
}

// Find the floor of n^{1/b}.
// In other words, the largest natural a such that a^b <= b
long floor_root(long n, long b) {
    long low = 2;
    long high = 1 << (((floor_log2(n) + 1) / b) + 1);
    long mid; long curr;
    while(1) {
        mid = (low + high) / 2;
        if (ipow(mid, b) > n)
            high = mid;
        else {
            low = mid;
            if(ipow(mid + 1, b) > n)
                return mid;
        }
    }
}

// Check if n is a^b with b > 2 and a natural by computing the floor of the bth
// root of n for all b up to log(n).

// TODO: Check if GMP's function for that (section 15.5.4 in the GMP PDF
// documentation) is at least as efficient. The docs say it's "not
// particularly efficient", but the below algorithm is not too sophisticated, so
// perhaps it's not as efficient as can be, but at least as efficient as this.
int is_perfect_power(long n) {
    for(long b = 2; b <= floor_log2(n); b++)
        if(ipow(floor_root(n, b), b) == n) return 1;
    return 0;
}

// Power n^k modulo r, TODO: replace with efficient GMP implementation
long modpow(long r, long n, long k) {
    if(n > r) return modpow(r, n % r, k);
    if(k == 0) return 1;
    if(k % 2 == 0) {
        long inter = modpow(r, n, k >> 1);
        return (inter * inter) % r;
    }
    return (modpow(r, n, k-1) * n) % r;
}

// a * b mod n assuming a, b <= n
// TODO: replace with GMP multiplication + mod
long modmul(long n, long a, long b) {
    return (a * b) % n;
}

// a + b mod n
long addmod(long n, long a, long b) {
    return (a + b) % n;
}


// TODO: Evaluate phi(r) and log n more rigourously for a tighter bound, if
// indeed found to be worthwhile to compute phi(r) (probably is, but needs
// testing).
long bound(long r, long n) {
    long rt = floor_root(r - 1, 2) + 1; // approximately ceil(sqrt(r - 1))
    return rt * (floor_log2(n) + 1);
}

// Find the minimum r for which n has (multiplicative) order > log^2 n in Z_r.
// Assumes: n is not a perfect power. So either n = 1, 2, or log(n) is not an
// integer, so floor_log2(n) < log(n). Return 0 if n is found to be composite at
// this step.
long min_r(long n) {
    long r;
    long lg = floor_log2(n);
    for(int r = 2; ;r++) {
        // TODO: replace with GMP gcd
        long d = gcd(r, n);
        if(1 < d && d < n)
            return 0;
        // Check if the order of n in Z_r is more than log^2.
        for(long k = 1; k <= log2(n) * log2(n); k++) { // TODO: take into account the
                                                       // accuracy of log2(n)
            if(modpow(r, n, k) == 1)
                goto next_r; // This r doesn't work, skip to the next one (outer
                             // continue)
        }
        // This r works, return it.
        return r;
        next_r:
    }
}
