#include <stdio.h>  // For: perror

int check(int n, double *A, double *B) {
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            if (A[i + j * n] != B[i + j * n]) {
                printf("%i %f, %f Diff  \n", n, A[i + j * n], B[i + j * n]);
                return 0;
            }
        }
    return 1;
}