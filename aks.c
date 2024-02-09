#include<stdio.h>
#include<stdlib.h>
#include<pthread.h>
#include<unistd.h>
#include<time.h>
#include"arithmetic.h"
#include"polynomial.h"
#define PROGRESS_INTERVAL 2
#define N_THREADS 4


// Check if the AKS equation for r, n, and a given.
int eqn(unsigned long r, const mpz_t n, const mpz_t a) {
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
    //print_poly(r, P);
    //print_poly(r, Q);
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
    mpz_t limit;
    mpz_init(limit);
    bound(limit, r, n);
    unsigned long limit_ui = mpz_get_ui(limit);
    time_t startTime = time(NULL);
    time_t currentTime;
    
    // Shared thread variable stating that there is still no proof that n is
    // composite.
    int still_prime = 1;
    pthread_t threads[N_THREADS];
    struct thread_arg {
        int thread_id;
        mpz_t a_mpz;
    };
    struct thread_arg thread_args[N_THREADS];

    void* thread_function(void* args) {
        struct thread_arg the_args = *((struct thread_arg*) args);
        int thread_id = the_args.thread_id;
        for(unsigned long i = 1; i <= limit_ui/N_THREADS; i++) { // AKS step 5
            // TODO: add the possibility of a --verbose option to control logging
            //       levels.
            unsigned long a = i * N_THREADS + thread_id + 1;
            currentTime = time(NULL);
            if(currentTime - startTime >= 2) {
                printf("\rProgress: %d/%d", a, limit_ui);
                fflush(stdout);
                startTime = currentTime;
            }
            if(!still_prime)
                return NULL;
            mpz_set_ui(the_args.a_mpz, a);
            if(!eqn(r, n, the_args.a_mpz)) {
                gmp_printf("thread %d: r, n, a = %d, %Zd, %d\n", thread_id, r, n, a);
                still_prime = 0;
                return NULL;
            }
        }
    }
    
    for(int i = 0; i < N_THREADS; i++) {
        thread_args[i].thread_id = i;
        mpz_init(thread_args[i].a_mpz);
        pthread_create(&threads[i], NULL, thread_function, &thread_args[i]);
    }

    for(int i = 0; i < N_THREADS; i++) {
        pthread_join(threads[i], NULL);
    }

    printf("\n");
    mpz_clear(limit);
    return still_prime;
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
    //printf("%d", is_prime(i_mpz));
    clock_t start, end;

    unsigned long r = 11;
    unsigned long n_P = 8;
    unsigned long n_Q = 3;
    mpz_set_ui(n, 29);
    mpz_t* P = (mpz_t*) malloc(r * sizeof(mpz_t));
    poly_init(r, P);
    mpz_t* Q = (mpz_t*) malloc(r * sizeof(mpz_t));
    poly_init(r, Q);
    mpz_set_ui(P[0], 4);
    mpz_set_ui(P[7], 1);
    mpz_set_ui(Q[0], 8);
    mpz_set_ui(Q[1], 2);
    mpz_set_ui(Q[2], 27);
    // P = 4 + X^7
    // Q = 8 + 2 X + 27 X^2
    // n = 29
    // r = 13
    // R = P * Q mod n mod X^r - 1 = 27*x^9 + 2*x^8 + 8*x^7 + 21*x^2 + 8*x + 3 (sage)
    // R_ = P^2 mod n mod X^r - 1 = 8*x^7 + x^3 + 16 (sage)
    mpz_t* R = (mpz_t*) malloc(r * sizeof(mpz_t));
    mpz_t* R_ = (mpz_t*) malloc(r * sizeof(mpz_t));
    poly_init(r, R);
    poly_init(r, R_);
    //shifted_add_mod_poly(r, n, R, n_P, P, n_Q, Q, 20);
    printf("P * Q starts here\n");
    karatsuba_poly_mul_mod(r, n, R, n_P, P, n_Q, Q, 0);
    karatsuba_poly_mul_mod(r, n, R_, n_P, P, n_P, P, 0);
    print_poly(r, P);
    print_poly(r, Q);
    print_poly(r, R);
    print_poly(r, R_);
    poly_free(r, P);
    poly_free(r, Q);
    poly_free(r, R);
    poly_free(r, R_);

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
        b = is_prime(i_mpz);
        end = clock();
        unsigned long r = min_r(i_mpz);
        mpz_t bnd;
        mpz_init(bnd);
        bound(bnd, r, n);
        printf("%ld r: %ld euler: %ld\n", mpz_get_ui(i_mpz), r, euler(r));
        //gmp_printf("%Zd prime: %d\tTime: %f r: \n", i_mpz, b, (float)(end - start)/CLOCKS_PER_SEC, min_r(i_mpz));
    }
        
    mpz_clear(i_mpz);
    mpz_clear(n);
    return 0;
}
