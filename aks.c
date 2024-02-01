#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<time.h>
#include"arithmetic.h"
#include"polynomial.h"


// Check if the AKS equation for r, n, and a given.
int eqn(unsigned long r, mpz_t n, const mpz_t a) {
    mpz_t a_; mpz_init(a_);
    mpz_mod(a_, a, n);
    mpz_t n_mod_r_mpz;
    mpz_init(n_mod_r_mpz);
    mpz_mod_ui(n_mod_r_mpz, n, r);
    unsigned long n_mod_r = mpz_get_ui(n_mod_r_mpz);
    // The following two lines are essentially
    // why r has to be an unsigned long, not mpz_t
    mpz_t* P = (mpz_t*) malloc(r * sizeof(mpz_t));
    mpz_t* Q = (mpz_t*) malloc(r * sizeof(mpz_t));
    poly_init(r, P); poly_init(r, Q);
    zero_poly(r, P);
    zero_poly(r, Q);
    mpz_set(P[0], a);
    mpz_set_ui(P[1], 1);
    poly_pow_mod(r, n, Q, 2, P, n);
    mpz_sub_ui(Q[n_mod_r], Q[n_mod_r], 1);
    mpz_sub(Q[0], Q[0], a);
    mpz_mod(Q[n_mod_r], Q[n_mod_r], n);
    mpz_mod(Q[0], Q[0], n);
    int b = is_zero_poly(r, Q);
    poly_free(r, P); poly_free(r, Q);
    mpz_clear(a_);
    return b;
}


// The actual AKS algorithm
int is_prime(mpz_t n) {
    if(mpz_perfect_power_p(n)) return 0; // AKS step 1
    unsigned long r = min_r(n); //AKS step 2
    if(r <= 1) return r; // AKS steps 3 and 4
    mpz_t limit, a_mpz, remainder;
    mpz_init(limit);
    mpz_init(a_mpz);
    mpz_init(remainder);
    bound(limit, r, n);
    unsigned long limit_ui = mpz_get_ui(limit);
    for(unsigned long a = 1; a < limit_ui; a++) { // AKS step 5
        // TODO: add the possibility of a --verbose option to control logging
        //       levels.
        mpz_set_ui(a_mpz, a);
        if(!eqn(r, n, a_mpz)) return 0;
    }
    mpz_clear(limit);
    mpz_clear(a_mpz);
    mpz_clear(remainder);
    //printf("\n")
    return 1;
}



int main() {
    mpz_t n;
    mpz_init(n);
    mpz_set_ui(n, 19997);
    mpz_t a;
    mpz_init(a);
    mpz_set_ui(a, 52);
    printf("%ld\n", min_r(n));
    //printf("%d\n", eqn(min_r(n), n, a));
    mpz_t i_mpz;
    mpz_init(i_mpz);
    int b;
    for(unsigned long int i = 1; i < 200; i++) {
        mpz_set_ui(i_mpz, i);
        //printf("n = %ld; r: %ld\n", i, min_r(i_mpz));
        b = is_prime(i_mpz);
        if(b)
            printf("n = %ld; prime: %d\n", i, b);
    }
    mpz_set_ui(i_mpz, 672629);
    mpz_set_ui(i_mpz, 542497301);
    mpz_set_ui(i_mpz, 54294967291);
    printf("%d", is_prime(i_mpz));
    clock_t start, end;

//    unsigned long known_primes[] = {2, 5, 11, 23, 47, 97, 197, 397, 797, 1597,
//        3203, 6421, 12853, 25717, 51437, 102877, 205759, 411527, 823117,
//        1646237, 3292489, 6584983, 13169977, 26339969, 52679969, 105359939,
//        210719881, 421439783, 842879579, 1685759167, 3371518343, 6743036717,
//        13486073473, 26972146961, 53944293929, 107888587883, 215777175787,
//        431554351609, 863108703229}; // 39

//    unsigned long known_primes[] = {3, 5, 11, 17, 37, 67, 131, 257, 521, 1031,
//        2053, 4099, 8209, 16411, 32771, 65537, 131101, 262147, 524309, 1048583,
//        2097169, 4194319, 8388617, 16777259, 33554467, 67108879, 134217757,
//        268435459, 536870923, 1073741827, 2147483659, 4294967311, 8589934609,
//        17179869209, 34359738421, 68719476767, 137438953481, 274877906951,
//        549755813911, 1099511627791, 2199023255579, 4398046511119}; // 42

    // On the edge of what unsigned long can still hold 
    //

    unsigned long known_primes[] = {3, 5, 11, 17, 37, 67, 131, 257, 521, 1031,
        2053, 4099, 8209, 16411, 32771, 65537, 131101, 262147, 524309, 1048583,
        2097169, 4194319, 8388617, 16777259, 33554467, 67108879, 134217757,
        268435459, 536870923, 1073741827, 2147483659, 4294967311, 8589934609,
        17179869209, 34359738421, 68719476767, 137438953481, 274877906951,
        549755813911, 1099511627791, 2199023255579, 4398046511119,
        8796093022237, 17592186044423, 35184372088891, 70368744177679,
        140737488355333, 281474976710677, 562949953421381, 1125899906842679}; // 50

    for(int i = 0; i < 50; i++) {
        mpz_set_ui(i_mpz, known_primes[i]);
        unsigned long r = min_r(i_mpz);
        printf("%ld \t r: %ld \t phi: %ld\n", known_primes[i], r, euler(r));
    }

    for(int i = 0; i < 50; i++) {
        mpz_set_ui(i_mpz, known_primes[i]);
        start = clock();
        b = is_prime(i_mpz);
        end = clock();
        printf("%ld %dms\n", known_primes[i], ((end - start)*1000)/CLOCKS_PER_SEC);
    }

    while(1) {
        printf("Input test number (base 10): ");
        gmp_scanf("%Zd", i_mpz);
        start = clock();
        //b = is_prime(i_mpz);
        end = clock();
        unsigned long r = min_r(i_mpz);
        mpz_t bnd;
        mpz_init(bnd);
        bound(bnd, r, n);
        printf("%ld r: %ld euler: %ld\n", mpz_get_ui(i_mpz), r, euler(r));
        //gmp_printf("%Zd prime: %d\tTime: %f r: \n", i_mpz, b, (float)(end - start)/CLOCKS_PER_SEC, min_r(i_mpz));
    }
        
    mpz_clear(i_mpz);
    return 0;
}
