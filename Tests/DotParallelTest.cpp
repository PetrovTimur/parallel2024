#include <cmath>
#include <omp.h>

#include "Solver/Kernels/mathfunc.h"


int main() {
    std::vector<double> x(1e8);
    std::vector<double> y(1e8);
    double res_single, res_multi;;

    for (std::size_t i = 0; i < x.size(); i++) {
        x[i] = std::sin(i);
        y[i] = std::cos(i);
    }

    const int max_threads = omp_get_max_threads();

    // Single-thread run
    omp_set_num_threads(1);
    dot(x.data(), y.data(), x.size(), res_single);

    //Multithread run
    omp_set_num_threads(max_threads);
    dot(x.data(), y.data(), x.size(), res_multi);

    return fabs(res_single - res_multi) > 1e-3;
}