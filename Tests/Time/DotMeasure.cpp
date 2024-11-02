#include <cmath>
#include <omp.h>
#include <iostream>
#include <ostream>
#include "Solver/Kernels/mathfunc.h"

int main() {
    omp_set_num_threads(omp_get_max_threads());

    std::vector<double> x;
    std::vector<double> y;
    double res;

    for (int k = 1e4; k <= 1e7; k *= 10) {
        x.resize(k);
        y.resize(k);

        int runs = 5;
        double aggregate_time = 0;
        for (int p = 0; p < runs; ++p) {
            // Fill
            #pragma omp parallel for
            for (int i = 0; i < k; i++) {
                x[i] = std::sin(i);
                y[i] = std::cos(i);
            }

            // Calculate
            double start = omp_get_wtime();
            dot(x, y, res);
            double end = omp_get_wtime();

            aggregate_time += end - start;

            std::cout << end - start << std::endl;
        }

        double average_time = aggregate_time / runs;

        std::cout << "N = " << k <<", average time: " << average_time << std::endl;
    }

}
