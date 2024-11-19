#include <cassert>
#include <cmath>
#include <omp.h>
#include <iostream>
#include <ostream>
#include <unistd.h>

#include "Solver/Kernels/mathfunc.h"

int main() {
    omp_set_num_threads(omp_get_max_threads());

    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> res;

    int T  = omp_get_max_threads();
    std::cout << "T = " << T << std::endl;

    for (int k = 1e5; k <= 1e8; k *= 10) {
        x.resize(k);
        y.resize(k);
        res.resize(k);

        // Fill
        #pragma omp parallel for proc_bind(spread)
        for (int i = 0; i < k; i++) {
            x[i] = sin(i);
            y[i] = cos(i);
        }

        int runs = 1e9 / k + 1;
        double aggregate_time = 0;
        for (int p = 1; p < runs + 1; ++p) {
            double alpha = std::tan(p);
            // Calculate
            double start = omp_get_wtime();
            axpy(alpha, x, y, res);
            double end = omp_get_wtime();

            aggregate_time += end - start;

        }

        double average_time = aggregate_time / runs;

        std::cout << 2 * k / (average_time * 1e9) << ", ";
    }
    std::cout << std::endl;
}
