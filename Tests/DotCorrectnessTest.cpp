#include <iostream>
#include <ostream>

#include "Solver/Kernels/mathfunc.h"


int main() {
    // Test 1
    std::vector<double> x = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
    std::vector<double> y = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};

    double res;
    dot(x, y, res);
    bool fl = res == 206;


    // Test 2
    x.resize(1e6);
    y.resize(1e6);

    for (int i = 0; i < x.size(); i++) {
        x[i] = 1. / (3 * (i + 1));
        y[i] = i + 1;
    }
    dot(x, y, res);

    fl = fl && (res == (x.size() / 3));

    return fl;
}