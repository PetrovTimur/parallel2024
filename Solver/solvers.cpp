#include "solvers.h"

#include <iostream>
#include <ostream>
#include "Kernels/mathfunc.h"

#ifdef USE_MPI
#include <mpi.h>
#endif

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
    double buf, total;
    int k = 0;

    do {
        k++;
        #pragma omp parallel for
        for (int i = 0; i < N; i++) {
            z[i] = r[i] / diag[i];
        }

        rho[0] = rho[1];

        dot(r, z, buf);
        #ifdef USE_MPI
        MPI_Allreduce(&buf, &rho[1], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        #else
        rho[1] = buf;
        #endif

        if (k == 1)
            #pragma omp parallel for
            for (int i = 0; i < N; i++)
                p[i] = z[i];
        else {
            double beta = rho[1] / rho[0];
            axpy(beta, p, z, p);
        }

        spMV(ia, ja, a, p, q);

        dot(p, q, buf);
        #ifdef USE_MPI
        MPI_Allreduce(&buf, &total, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        double alpha = rho[1] / total;
        #else
        double alpha = rho[1] / buf;
        #endif

        axpy(alpha, p, x, x);
        axpy(-alpha, q, r, r);
    }
    while (rho[1] > eps * eps && k < maxit);

    #pragma omp parallel for
    for (int i = 0; i < N; i++)
        res[i] = x[i];

    // std::cout << "k = " << k << std::endl;
    return k;
}
