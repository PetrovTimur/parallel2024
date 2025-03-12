#include <cmath>
#include <omp.h>

#include "Solver/Kernels/mathfunc.h"


int main() {
    std::vector<double> x(1e8);
    std::vector<double> y(1e8);
    std::vector<double> z(1e8);
    double norm_single, norm_multi;
    double alpha = 0.5163;

    for (std::size_t i = 0; i < x.size(); i++) {
        x[i] = std::sin(i);
        y[i] = std::cos(i);
    }

    int max_threads = omp_get_max_threads();

    // Single-thread run
    omp_set_num_threads(1);
    axpy(alpha, x.data(), y.data(), x.size(), z.data());
    dot(z.data(), z.data(), z.size(), norm_single);

    //Multithread run
    omp_set_num_threads(max_threads);
    axpy(alpha, x.data(), y.data(), x.size(), z.data());
    dot(z.data(), z.data(), z.size(), norm_multi);

    return fabs(norm_single - norm_multi) > 1e-3;
}