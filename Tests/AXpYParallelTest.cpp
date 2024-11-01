#include <cmath>
#include <omp.h>

#include "Solver/Kernels/mathfunc.h"


int main() {
    std::vector<double> x(1e8);
    std::vector<double> y(1e8);
    std::vector<double> z(1e8);
    double norm_single, norm_multi;
    double alpha = 0.5163;

    for (int i = 0; i < x.size(); i++) {
        x[i] = std::sin(i);
        y[i] = std::cos(i);
    }

    // Single-thread run
    omp_set_num_threads(1);
    axpy(alpha, x, y, z);
    dot(z, z, norm_single);

    //Multithread run
    omp_set_num_threads(omp_get_max_threads());
    axpy(alpha, x, y, z);
    dot(z, z, norm_multi);

    return norm_single != norm_multi;
}