#include "scan.h"
#include <omp.h>
#include <vector>

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
