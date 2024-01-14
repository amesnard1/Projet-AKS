#include<stdio.h>
#include<stdlib.h>
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

// To be replaced with GMP functions
long ipow_rec(long a, long b) {
    if(b == 0) return 1;
    if(b % 2 == 0) {
        long inter = ipow_rec(a, b >> 1);
        return inter * inter;
    }
    return ipow_rec(a, b-1) * a;
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

// TODO: Evaluate phi(r) and log n more rigourously for a tighter bound, if
// indeed found to be worthwhile to compute phi(r) (probably is, but needs
// testing).
long bound(long r, long n) {
    long rt = floor_root(r - 1, 2) + 1; // approximately ceil(sqrt(r - 1))
    return rt * (floor_log2(n) + 1);
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

// Check if the AKS equation for r, n, and a given.
int eqn(long r, long n, long a) {
    a = a % n;
    long* P = (long*) malloc(r * sizeof(long));
    long* Q = (long*) malloc(r * sizeof(long));
    for(long i = 0; i < r; i++) {
        P[i] = 0;
        Q[i] = 0;
    }
    P[0] = a;
    P[1] = 1;
    poly_pow_mod(r, n, Q, P, n);
    Q[n % r] -= 1;
    Q[0] -= a;
    Q[n % r] = Q[n % r] % n;
    Q[0] = Q[0] % n;
    int b = is_zero_poly(r, Q);
    free(P); free(Q);
    return b;
}


// The actual AKS algorithm
int is_prime(long n) {
    if(is_perfect_power(n)) return 0; // AKS step 1
    long r = min_r(n); // AKS step 2 (and most of 3)
    // For the next line, recall that min_r returns 0 if it detects a proper
    // factor of n before it finds an r such that o_r(n) > log^2 n.
    if(r == 0) return 0; // AKS step 3
    if(n <= r) return 1; // AKS step 4
    for(long a = 1; a <= bound(r, n); a++) // AKS step 5
        if(!eqn(r, n, a)) return 0;
    return 1;
}



int main() {
    printf("%d\n", floor_log2(20));
    printf("%d\n", ipow(2, 4));
    printf("%d\n", ipow(3, 5));
    printf("%d\n", floor_root(1241, 5));
    printf("%d\n", is_perfect_power(1419856));
    printf("%d\n", is_perfect_power(1419857));
    printf("%d\n", min_r(19997));
    printf("%d\n", eqn(min_r(19997), 19997, 52));
    long r = 10;
    long* P = (long*) malloc(r * sizeof(long));
    long* Q = (long*) malloc(r * sizeof(long));
    long* R = (long*) malloc(r * sizeof(long));
    zero_poly(r, P);
    zero_poly(r, Q);
    zero_poly(r, R);
    P[0] = 7; P[1] = 5; P[2] = 3; P[3] = 2;
    Q[0] = 1; Q[1] = 8; Q[2] = 4;
    poly_mul_mod(r, 100, R, P, Q);
    print_poly(r, R, r-1);
    zero_poly(r, R);
    poly_pow_mod(r, 18, R, P, 15);
    print_poly(r, R, r-1);
    for(int i = 1; i < 200; i++)
        if(is_prime(i))
            printf("n = %d; prime: %d\n", i, is_prime(i));
    return 0;
}
