#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"arithmetic.h"
#include"polynomial.h"


// Check if the AKS equation for r, n, and a given.
int eqn(long r, long n, long a) {
    a = a % n;
    long* P = (long*) malloc(r * sizeof(long));
    long* Q = (long*) malloc(r * sizeof(long));
    zero_poly(r, P);
    zero_poly(r, Q);
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
    // First eliminate trivial cases that do not fit into the algorithm
    if(n == 0 || n == 1) return 0;
    if(is_perfect_power(n)) return 0; // AKS step 1
    long r = min_r(n); // AKS step 2 (and most of 3)
    // For the next line, recall that min_r returns 0 if it detects a proper
    // factor of n before it finds an r such that o_r(n) > log^2 n.
    if(r == 0) return 0; // AKS step 3
    if(n <= r) return 1; // AKS step 4
    for(long a = 1; a <= bound(r, n); a++) { // AKS step 5
// Later, the following will be added conditional upon a --verbose flag from
// shell.
//        if(a % 10 == 0)
//            printf("progress: %d/%d\n", a, bound(r, n));
        if(!eqn(r, n, a)) return 0;
    }
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
    printf("%d", is_prime(672629));
    return 0;
}
