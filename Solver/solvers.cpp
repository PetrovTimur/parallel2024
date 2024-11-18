#include "solvers.h"

#include <iostream>
#include <ostream>

#include "Kernels/mathfunc.h"

int solve(std::vector<int> &ia, std::vector<int> &ja, std::vector<double> &a, std::vector<double> &b,
          std::vector<double> &diag, std::vector<double> &res) {
    int N = ia.size() - 1;
    double eps = 1e-3;
    int maxit = 100000;

    std::vector<double> x(N);
    std::vector<double> r(b);
    std::vector<double> z(N);
    std::vector<double> p(N);
    std::vector<double> q(N);
    std::vector<double> rho(2);
    double buf;
    int k = 0;

    do {
        k++;
        #pragma omp parallel for proc_bind(spread)
        for (int i = 0; i < N; i++) {
            z[i] = r[i] / diag[i];
        }

        rho[0] = rho[1];
        dot(r, z, rho[1]);

        if (k == 1)
            #pragma omp parallel for proc_bind(spread)
            for (int i = 0; i < N; i++)
                p[i] = z[i];
        else {
            double beta = rho[1] / rho[0];
            axpy(beta, p, z, p);
        }

        spMV(ia, ja, a, p, q);
        dot(p, q, buf);
        double alpha = rho[1] / buf;

        axpy(alpha, p, x, x);
        axpy(-alpha, q, r, r);
    }
    while (rho[1] > eps * eps && k < maxit);

    #pragma omp parallel for proc_bind(spread)
    for (int i = 0; i < N; i++)
        res[i] = x[i];

    // std::cout << "k = " << k << std::endl;
    return k;
}
