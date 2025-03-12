#include "Solver/Kernels/mathfunc.h"


int main() {
    // Test 1
    std::vector<double> x = {2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0};
    std::vector<double> y = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
    std::vector<double> res(x.size());
    double alpha = -0.5, norm;

    axpy(alpha, x.data(), y.data(), x.size(), res.data());
    dot(res.data(), res.data(), res.size(), norm);
    bool fl = norm == 0;


    // Test 2
    x.resize(1e6);
    y.resize(1e6);
    res.resize(1e6);
    alpha = 373;

    for (std::size_t i = 0; i < x.size(); i++) {
        x[i] = 1. / (alpha * (i + 1));
        y[i] = i + 1;
    }
    axpy(-alpha, x.data(), y.data(), x.size(), res.data());
    dot(res.data(), res.data(), res.size(), norm);
    fl = fl && norm == 0;

    return fl;
}