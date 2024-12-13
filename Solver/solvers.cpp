#include "solvers.h"

#include <iostream>
#include <ostream>
#include <unistd.h>

#include "Kernels/mathfunc.h"
#ifdef USE_MPI
#include "Utilities/coms.h"
#endif

#ifdef USE_MPI
#include <mpi.h>
#endif

#ifdef USE_MPI
int solve(int MyID, int Px, int top_halo, int left_halo, int right_halo, int bottom_halo, int i_count, int j_count,
        std::vector<int> &ia, std::vector<int> &ja, std::vector<double> &a, std::vector<double> &b,
        std::vector<double> &diag, std::vector<double> &res) {

    std::vector<double> recv_buf;
    std::vector<int> recv_offset(7);

    std::vector<double> send_buf;
    std::vector<int> send_offset(7);

    ComInitOffsets(top_halo, left_halo, right_halo, bottom_halo, i_count, j_count, recv_offset, send_offset);
    recv_buf.resize(recv_offset[recv_offset.size() - 1]);
    send_buf.resize(send_offset[send_offset.size() - 1]);

    int N = ia.size() - 1;
    double eps = 1e-3;
    int maxit = 20;

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
        #pragma omp parallel for proc_bind(master)
        for (int i = 0; i < N; i++) {
            z[i] = r[i] / diag[i];
        }

        rho[0] = rho[1];

        dot(r, z, buf);
        MPI_Allreduce(&buf, &rho[1], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        if (k == 1)
            #pragma omp parallel for proc_bind(master)
            for (int i = 0; i < N; i++)
                p[i] = z[i];
        else {
            double beta = rho[1] / rho[0];
            axpy(beta, p, z, p);
        }

        Com(MyID, Px, top_halo, left_halo, right_halo, bottom_halo, i_count, j_count,
        p, recv_offset, send_offset,
        recv_buf, send_buf);

        spMV(ia, ja, a, p, recv_buf, q);

        dot(p, q, buf);
        MPI_Allreduce(&buf, &total, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        double alpha = rho[1] / total;

        axpy(alpha, p, x, x);
        axpy(-alpha, q, r, r);
    }
    while (rho[1] > eps * eps && k < maxit);

    #pragma omp parallel for proc_bind(master)
    for (int i = 0; i < N; i++)
        res[i] = x[i];

    // std::cout << "k = " << k << std::endl;
    return k;
}
#else
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

        #ifdef USE_DEBUG_MODE
        std::cout << " k = " << k << std::endl;
        std::cout << " rho = " << rho[0] << ", " << rho[1] << std::endl;
        std::cout << " alpha = " << alpha << std::endl;
        #endif

    }
    while (rho[1] > eps * eps && k < maxit);

    #pragma omp parallel for proc_bind(spread)
    for (int i = 0; i < N; i++)
        res[i] = x[i];

    // std::cout << "k = " << k << std::endl;
    return k;
}
#endif