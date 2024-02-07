#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>

#define MAX_PRIMES 10000

//Function to generate the first 10,000 prime numbers using the Sieve of Eratosthenes (crible d'eratosthenes)

void generate_primes(int primes[], int max_primes) {
    bool is_prime[max_primes + 1];
    for (int i = 0; i <= max_primes; ++i) {
        is_prime[i] = true;
    }

    int count = 0;
    for (int p = 2; p <= max_primes; ++p) {
        if (is_prime[p]) {
            primes[count++] = p;
            for (int i = p * 2; i <= max_primes; i += p) {
                is_prime[i] = false;
            }
        }
    }
}

//function to test our function on different primes and give the time associated to each value
int main() {
    int prime_numbers[MAX_PRIMES];
    generate_primes(prime_numbers, MAX_PRIMES);
    for (int i = 0; i < MAX_PRIMES; ++i) {
        clock_t start_time, end_time;
        start_time = clock();
        ALGOAKS(prime_numbers[i])
        end_time = clock();
        printf("%f\n", (start_time - end_time))
    }

    return 0;
}
