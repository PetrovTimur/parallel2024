#include "solvers.h"

#include <algorithm>
#include <complex>
#include <iostream>
#include <omp.h>
#include <ostream>
#include <unistd.h>
#include <Utilities/logger.h>

#include "Kernels/mathfunc.h"
#ifdef USE_MPI
#include "Utilities/coms.h"
#include <mpi.h>
#endif

#ifdef USE_MPI
int solve(const int MyID, const std::vector<int> &Part, const std::vector<int> &L2G, const std::vector<int> &ia, const std::vector<int> &ja,
          const std::vector<double> &a, const std::vector<double> &b, const std::vector<double> &diag, std::vector<double> &res, const double eps, int maxit) {

    std::vector<int> RecvOffset, SendOffset;
    std::vector<int> Recv, Send;
    std::vector<int> Neighbors;

    ComInit(ia, ja, Part, L2G, MyID, RecvOffset, SendOffset, Send, Recv, Neighbors);

    // std::cout << "My ID: " << MyID << ", Neighbors size: " << Neighbors.size() << std::endl;

    int N = ia.size() - 1;
    // const double eps = 1e-3;
    // const int maxit = 20;

    const auto z = new double[N];
    const auto p = new double[L2G.size()];

    const auto q = new double[N];
    arrInit(q, 0., N);

    const auto x = new double[N];
    arrInit(x, 0., N);

    const auto r = new double[N];
    arrCopy(r, b.data(), N);

    const auto rho = new double[2];
    rho[0] = rho[1] = 0;

    double buf, total, norm;
    int k = 0;

    do {
        k++;

        if (MyID == 0) {
            LOG_INFO << "Iteration " << k << std::endl;
        }

        MPI_Barrier(MPI_COMM_WORLD);
        double start = MPI_Wtime();

        // #pragma omp parallel for proc_bind(master)
        for (int i = 0; i < N; i++) {
            z[i] = r[i] * diag[i];
        }

        rho[0] = rho[1];

        dot(r, z, N, buf);
        MPI_Allreduce(&buf, &rho[1], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        if (k == 1)
            // #pragma omp parallel for proc_bind(master)
            for (int i = 0; i < N; i++)
                p[i] = z[i];
        else {
            double beta = rho[1] / rho[0];
            axpy(beta, p, z, N, p);
        }

        ComUpdate(p, RecvOffset, SendOffset, Neighbors, Send, Recv);

        spMV(ia.data(), ja.data(), a.data(), p, ia.size(), q);

        dot(p, q, N, buf);
        MPI_Allreduce(&buf, &total, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        double alpha = rho[1] / total;

        axpy(alpha, p, x, N, x);
        axpy(-alpha, q, r, N, r);

        MPI_Barrier(MPI_COMM_WORLD);
        double end = MPI_Wtime();
        dot(x, x, N, norm);
        MPI_Allreduce(&norm, &total, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        if (MyID == 0) {
            LOG_INFO << "Time " << end - start << std::endl;
            // LOG_INFO << "rho = " << rho[0] << ", " << rho[1] << ", alpha = " << alpha << std::endl;
            LOG_INFO << "L2 norm: " << std::sqrt(total) << std::endl;
            LOG << "--------------------------------------------------" << std::endl << std::endl;
        }
    }
    while (rho[1] > eps * eps && k < maxit);

    // #pragma omp parallel for proc_bind(master)
    for (int i = 0; i < N; i++)
        res[i] = x[i];

    delete[] z;
    delete[] p;
    delete[] q;
    delete[] x;
    delete[] r;
    delete[] rho;
    return k;
}
#else
int solve(const int *ia, const int *ja, const double *a, const double *b, const double *diag, const int size, double *res, const double eps, const int maxIterations) {
    const int N = size;

    const auto z = new double[N];
    const auto p = new double[N];

    const auto q = new double[N];
    arrInit(q, 0., N);

    const auto x = new double[N];
    arrInit(x, 0., N);

    const auto r = new double[N];
    arrCopy(r, b, N);

    const auto rho = new double[2];
    rho[0] = rho[1] = 0;

    double buf, norm;
    int k = 0;

    do {
        k++;

        LOG_INFO << "Iteration " << k << std::endl;

        const auto start = omp_get_wtime();

        #pragma omp parallel for default(none) shared(z, r, diag, N)
        for (int i = 0; i < N; i++) {
            z[i] = r[i] * diag[i];
        }

        rho[0] = rho[1];

        dot(r, z, N, rho[1]);

        if (k == 1)
            #pragma omp parallel for default(none) shared(p, z, N)
            for (int i = 0; i < N; i++)
                p[i] = z[i];
        else {
            const double beta = rho[1] / rho[0];
            axpy(beta, p, z, N, p);
        }

        spMV(ia, ja, a, p, N + 1, q);

        dot(p, q, N, buf);
        double alpha = rho[1] / buf;

        axpy(alpha, p, x, N, x);
        axpy(-alpha, q, r, N, r);

        LOG_INFO << "Time " << omp_get_wtime() - start << std::endl;
        // LOG_INFO << "rho = " << rho[0] << ", " << rho[1] << ", alpha = " << alpha << std::endl;
        dot(x, x, N, norm);
        LOG_INFO << "Solution norm: " << std::sqrt(norm) << std::endl;
        #ifdef DEBUG_MODE
        spMV(ia, ja, a, x, N + 1, z);

        #pragma omp parallel for default(none) shared(z, b, N)
        for (int i = 0; i < N; i++)
            z[i] -= b[i];

        dot(z, z, N, norm);
        LOG_DEBUG << "Residual norm: " << std::sqrt(norm) << std::endl;
        #endif

        LOG << "--------------------------------------------------" << std::endl << std::endl;

    }
    while (rho[1] > eps * eps && k < maxIterations);

    #pragma omp parallel for default(none) shared(res, x, N)
    for (int i = 0; i < N; i++)
        res[i] = x[i];

    delete[] z;
    delete[] p;
    delete[] q;
    delete[] x;
    delete[] r;
    delete[] rho;

    return k;
}
#endif