#include "solvers.h"
#include "Ops/mathfunc.h"

void solve(std::vector<int> &ia, std::vector<int> &ja, std::vector<double> &a, std::vector<double> &b,
    std::vector<double> &diag, std::vector<double> &res) {
    int N = ia.size() - 1;
    double eps = 0.001;
    int maxit = 100000;

    std::vector<double> x(N);
    std::vector<double> r(b);
    std::vector<double> z(N);
    std::vector<double> p(N);
    std::vector<double> q(N);
    std::vector<double> rho(2);
    int k = 0;

    do {
        k++;
        #pragma omp parallel for
        for (int i = 0; i < N; i++) {
            z[i] = r[i] / diag[i];
        }

        rho[0] = rho[1];
        rho[1] = dot(r, z);

        if (k == 1)
            #pragma omp parallel for
            for (int i = 0; i < N; i++)
                p[i] = z[i];
        else {
            double beta = rho[1] / rho[0];
            #pragma omp parallel for
            for (int i = 0; i < N; i++)
                p[i] = z[i] + beta * p[i];
        }

        spMV(ia, ja, a, p, q);
        double alpha = rho[1] / dot(p, q);

        #pragma omp parallel for
        for (int i = 0; i < N; i++)
            x[i] = x[i] + alpha * p[i];

        #pragma omp parallel for
        for (int i = 0; i < N; i++)
            r[i] = r[i] - alpha * q[i];
    }
    while (rho[1] > eps * eps && k < maxit);

    #pragma omp parallel for
    for (int i = 0; i < N; i++)
        res[i] = x[i];
}
