#include "mathfunc.h"

#include <iostream>
#include <omp.h>


void dot(const double *x, const double *y, const int size, double &res) {
    double sum = 0;

    #pragma omp parallel for default(none) shared(x, y, size) reduction(+:sum)
    for (int i = 0; i < size; i++) {
        sum += x[i] * y[i];
    }
    res = sum;
}

void axpy(const double a, const double *x, const double *y, const int size, double *res) {
    #pragma omp parallel for default(none) shared(a, x, y, size, res)
    for (int i = 0; i < size; i++) {
        res[i] = a * x[i] + y[i];
    }
}

void spMV(const int *ia, const int *ja, const double *a, const double *b, const int size, double *res) {
    #pragma omp parallel for default(none) shared(ia, ja, a, b, size, res)
    for (int i = 0; i < size - 1; i++) {
        double sum = 0;
        for (int col = ia[i]; col < ia[i + 1]; col++) {
            const int j = ja[col];
            const double a_ij = a[col];
            sum += a_ij * b[j];
        }
        res[i] = sum;
    }
}

#include <vector>
#include <algorithm>

void scan(const int* input, int* output, const int n) {
    int num_threads;
    #pragma omp parallel default(none) shared(num_threads)
    {
        #pragma omp single
        {
            num_threads = omp_get_num_threads();
        }
    }

    std::vector<int> partial(num_threads + 1, 0);

    #pragma omp parallel default(none) shared(partial, num_threads, input, output, n)
    {
        int tid = omp_get_thread_num();
        int chunk = (n + num_threads - 1) / num_threads;
        int start = tid * chunk;
        int end   = std::min(start + chunk, n);
        int sum = 0;
        for (int i = start; i < end; ++i) {
            sum      += input[i];
            output[i] = sum;
        }
        partial[tid + 1] = sum;

        #pragma omp barrier

        #pragma omp single
        {
            for (int i = 1; i <= num_threads; ++i)
                partial[i] += partial[i - 1];
        }

        #pragma omp barrier

        int offset = partial[tid];
        for (int i = start; i < end; ++i) {
            output[i] += offset;
        }
    }
}