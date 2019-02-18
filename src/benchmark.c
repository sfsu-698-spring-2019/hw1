#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#ifdef GETTIMEOFDAY
#include <sys/time.h>
#else

#include <time.h>

#endif

extern void square_dgemm(int, double *, double *, double *);

extern int check(int, double *, double *);

extern void square_dgemm_naive(int, double *, double *, double *);

double wall_time() {
#ifdef GETTIMEOFDAY
    struct timeval t;
    gettimeofday (&t, NULL);
    return 1.*t.tv_sec + 1.e-6*t.tv_usec;
#else
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return 1. * t.tv_sec + 1.e-9 * t.tv_nsec;
#endif
}

void die(const char *message) {
    perror(message);
    exit(EXIT_FAILURE);
}

void fill(double *p, int n) {
    for (int i = 0; i < n; ++i)
        p[i] = 2 * drand48() - 1; // Uniformly distributed over [-1, 1]
}

void copy(double *from, double *to, int n) {
    for (int i = 0; i < n; i++) {
        to[i] = from[i];
    }
}

int main() {
    int test_sizes[] = {31, 32, 96, 97, 127, 128, 129, 191, 192, 229, 255, 256, 257,
                        319, 320, 321, 417, 479, 480, 511, 512, 639, 640, 767, 768, 769};

    int nsizes = sizeof(test_sizes) / sizeof(test_sizes[0]);

    /* assume last size is also the largest size */
    int nmax = test_sizes[nsizes - 1];

    double *buf = NULL;
    double *variation_1 = (double *) malloc(nmax * nmax * sizeof(double));
    double *variation_2 = (double *) malloc(nmax * nmax * sizeof(double));
    double *variation_3 = (double *) malloc(nmax * nmax * sizeof(double));
    double *variation_4 = (double *) malloc(nmax * nmax * sizeof(double));
    double *variation_5 = (double *) malloc(nmax * nmax * sizeof(double));

    buf = (double *) malloc(3 * nmax * nmax * sizeof(double));
    if (buf == NULL) die("failed to allocate largest problem size");

    /* For each test size */
    double start_time, time, naive_time;
    int iterations = 10;
    for (int isize = 0; isize < nsizes; ++isize) {
        int n = test_sizes[isize];
        double *A = buf + 0;
        double *B = A + nmax * nmax;
        double *C = B + nmax * nmax;
        int square_n = n * n;

        fill(A, square_n);
        fill(B, square_n);
        fill(C, square_n);

        copy(C, variation_1, square_n);
        copy(C, variation_2, square_n);
        copy(C, variation_3, square_n);
        copy(C, variation_4, square_n);
        copy(C, variation_5, square_n);
        check(n, C, variation_1);

        for (int i = 0; i < iterations; i++) {
            start_time = wall_time();
            square_dgemm_naive(n, A, B, C);
            naive_time = wall_time() - start_time;
        }
        naive_time = naive_time / iterations;

        // Make one of these per variation
        for (int i = 0; i < iterations; i++) {
            start_time = wall_time();
            square_dgemm(n, A, B, variation_1);
            time = wall_time() - start_time;
        }
        // make sure the output is correct, but don't measure time of this
        if (check(n, C, variation_1) == 0) {
            printf("INVALID OUTPUT");
            exit(1);
        }
        time = time / iterations;

        printf("Variation 1: Naive: %f, Blocked: %f, Size: %i  \n", naive_time, time, n);
    }

    free(buf);
    free(variation_1);
    free(variation_2);
    free(variation_3);
    free(variation_4);
    free(variation_5);
    return 0;
}