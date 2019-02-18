#include "globals.h"
#include <stdio.h>

#define min(a, b) (((a)<(b))?(a):(b))

/* This auxiliary subroutine performs a smaller dgemm operation
 *  C := C + A * B
 * where C is M-by-N, A is M-by-K, and B is K-by-N. */
static void do_block(int lda, int M, int N, int K, double *A, double *B, double *C) {
    /* For each row i of A */
    for (int i = 0; i < M; ++i)
        /* For each column j of B */
        for (int j = 0; j < N; ++j) {
            /* Compute C(i,j) */
            double cij = C[i + j * lda];
            for (int k = 0; k < K; ++k)
                cij += A[i + k * lda] * B[k + j * lda];
            C[i + j * lda] = cij;
        }
}

static void do_block_v1(int lda, int M, int N, int K, double *A, double *B, double *C) {
    // hardcoded to only work with blocks specifically of size BLOCK_SIZE x BLOCK_SIZE
    static double a[BLOCK_SIZE * BLOCK_SIZE];

    // Make a cached copy of A so that each leap will now be k * BLOCK_SIZE
    // However it takes overhead to make the copy, so can be good or bad depending on your cache sizes
    for (int j = 0; j < K; j++)
        for (int i = 0; i < M; i++)
            a[i + j * BLOCK_SIZE] = A[i + j * lda];

    /* For each row i of A */
    for (int i = 0; i < M; ++i)
        /* For each column j of B */
        for (int j = 0; j < N; ++j) {
            /* Compute C(i,j) */
            double cij = C[i + j * lda];
            for (int k = 0; k < K; ++k)
                cij += a[i + k * BLOCK_SIZE] * B[k + j * lda];
            C[i + j * lda] = cij;
        }
}

/* This routine performs a dgemm operation
 *  C := C + A * B
 * where A, B, and C are lda-by-lda matrices stored in column-major format. 
 * On exit, A and B maintain their input values. */
void square_dgemm(int lda, double *A, double *B, double *C) {
    /* For each block-row of A */
    for (int i = 0; i < lda; i += BLOCK_SIZE)
        /* For each block-column of B */
        for (int j = 0; j < lda; j += BLOCK_SIZE)
            /* Accumulate block dgemms into block of C */
            for (int k = 0; k < lda; k += BLOCK_SIZE) {
                /* Correct block dimensions if block "goes off edge of" the matrix */
                int M = min (BLOCK_SIZE, lda - i);
                int N = min (BLOCK_SIZE, lda - j);
                int K = min (BLOCK_SIZE, lda - k);
                // can this be improved even more?
                if ((M % BLOCK_SIZE == 0) &&
                    (N % BLOCK_SIZE == 0) &&
                    (K % BLOCK_SIZE == 0)) {
                    do_block_v1(lda, M, N, K, A + i + k * lda, B + k + j * lda, C + i + j * lda);
                } else {
                    // this is the stock version, can you make this better?
                    do_block(lda, M, N, K, A + i + k * lda, B + k + j * lda, C + i + j * lda);
                }
            }
}