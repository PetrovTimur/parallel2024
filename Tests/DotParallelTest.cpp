#include <cmath>
#include <omp.h>

#include "Solver/Kernels/mathfunc.h"


int main() {
    std::vector<double> x(1e8);
    std::vector<double> y(1e8);
    double res_single, res_multi;;

    for (int i = 0; i < x.size(); i++) {
        x[i] = std::sin(i);
        y[i] = std::cos(i);
    }

    // Single-thread run
    omp_set_num_threads(1);
    dot(x, y, res_single);

    //Multithread run
    omp_set_num_threads(omp_get_max_threads());
    dot(x, y, res_multi);

    return res_single != res_multi;
}