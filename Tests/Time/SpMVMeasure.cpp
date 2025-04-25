#include <cassert>
#include <omp.h>
#include <iostream>
#include <ostream>
#include <tuple>

#include "Solver/csr.h"
#include "Solver/Kernels/mathfunc.h"
#include "Utilities/input.h"

#ifdef USE_MPI
#include <mpi.h>
#include "Utilities/coms.h"
#endif

int main(int argc, char *argv[]) {
    #ifdef USE_MPI
    omp_set_num_threads(omp_get_max_threads());

    int mpi_res;
    mpi_res = MPI_Init(&argc, &argv); // первым делом подключаемся к MPI
    if(mpi_res != MPI_SUCCESS) {
        MPI_Abort(MPI_COMM_WORLD, -1);
        exit(mpi_res);
    }

    int NumProc;
    mpi_res = MPI_Comm_size(MPI_COMM_WORLD,&NumProc); // узнаем число процессов
    if(mpi_res!= MPI_SUCCESS){
        MPI_Abort(MPI_COMM_WORLD, -1);
        exit(mpi_res);
    }

    int MyID;
    mpi_res = MPI_Comm_rank(MPI_COMM_WORLD,&MyID); // узнаем номер данного процесса
    if(mpi_res!= MPI_SUCCESS){
        MPI_Abort(MPI_COMM_WORLD, -1);
        exit(mpi_res);
    }

    if (MyID == 0) {
        std::cout << "NumProc = " << NumProc << std::endl;
        int T  = omp_get_max_threads();
        std::cout << "T = " << T << std::endl;
    }

    for (int k = 1e6; k <= 1e8; k *= 10) {
        int K1 = 30;
        int K2 = 23;
        int Nx = 1000;
        int Ny = k / Nx / 5;
        int Px = 2;
        int Py = NumProc / 2;

        std::vector<int> L2G;
        std::vector<int> Part;

        input(Nx, Ny, Px, Py, MyID, L2G, Part);

        std::vector<int> ia_en;
        std::vector<int> ja_en;

        makeIncidenceMatrixCSR(Nx, Ny, K1, K2, L2G, ia_en, ja_en, Part);
        std::unordered_map<int, int> G2L;
        fillG2L(L2G, G2L);

        std::unordered_map<int, int> G2L_nodes;
        constructG2L(ia_en, ja_en, G2L_nodes);

        localizeCSR(ia_en.data(), ia_en.size(), ja_en.data(), G2L_nodes);

        int *ia_ne, *ja_ne;
        transposeCSR(ia_en, ja_en, G2L_nodes.size(), ia_ne, ja_ne);

        std::vector<int> ia_ee, ja_ee;
        buildAdjacencyMatrixCSRUsingSort(ia_en.data(), ja_en.data(), ia_ne, ja_ne, ia_ee, ja_ee, ia_en.size() - 1, Part, MyID);

        std::vector<double> a(ja_ee.size());
        std::vector<double> b(ia_ee.size() - 1);
        std::vector<double> diag(ia_ee.size() - 1);

        fillCSR(ia_ee, ja_ee, L2G, a, b, diag);

        std::vector<double> res(ia_ee.size() - 1);

        // std::vector<int> RecvOffset, SendOffset;
        // std::vector<int> Recv, Send;
        // std::vector<int> Neighbors;

        // std::cout << "Ready " << ia_ee.size() - 1 << std::endl;

        int runs = 1e7 / k + 1;
        double aggregate_time = 0;
        for (int p = 0; p < runs; ++p) {
            // Calculate
            // ComInit(ia_ee, ja_ee, Part, L2G, MyID, RecvOffset, SendOffset, Send, Recv, Neighbors);
            std::vector<double> q(b.size() * 20, 1);
            double start = MPI_Wtime();
            spMV(ia_ee.data(), ja_ee.data(), a.data(), q.data(), ia_ee.size(), res.data());
            MPI_Barrier(MPI_COMM_WORLD);
            double end = MPI_Wtime();

            aggregate_time += end - start;
        }

        double average_time = aggregate_time / runs;

        if (MyID == 0)
            std::cout << 2 * NumProc * ia_ee.back() / (average_time * 1e9) << ",";
    }
    if (MyID == 0)
        std::cout << std::endl;

    MPI_Finalize();

    #else
    omp_set_num_threads(omp_get_max_threads());
    // std::vector<int> ia;
    // std::vector<int> ja;
    // std::vector<double> a;
    // std::vector<double> b;
    // std::vector<double> diag;
    // std::vector<double> res;

    int T  = omp_get_max_threads();
    std::cout << "T = " << T << std::endl;

    for (int k = 1e6; k <= 1e8; k *= 10) {
        int K1 = 30;
        int K2 = 23;
        int Nx = 1000;
        int Ny = k / Nx / 5;

        const gridInfo grid(Nx, Ny, K1, K2);
        int Nn = grid.totalNodes;
        int Ne = grid.totalElements;

        std::vector<int> ia_en;
        std::vector<int> ja_en;
        makeIncidenceMatrixCSR(Nx, Ny, K1, K2, 0, Ny - 1, 0, Nx - 1, ia_en, ja_en);

        int *ia_ne, *ja_ne;
        transposeCSR(ia_en, ja_en, Nn, ia_ne, ja_ne);

        auto matrix_ee = buildAdjacencyMatrixCSRUsingSort(ia_en.data(), ja_en.data(), ia_ne, ja_ne, Ne, Nn);
        int *ia = std::get<0>(matrix_ee);
        int *ja = std::get<1>(matrix_ee);

        #ifdef USE_CUDA
        auto a = new float[ia[Ne]];
        auto b = new float[Ne];
        auto diag = new float[Ne];
        #else
        auto a = new double[ia[Ne]];
        auto b = new double[Ne];
        auto diag = new double[Ne];
        #endif

        fillCSR(ia, ja, a, b, diag, Ne);

        #ifdef USE_CUDA
        auto res = new float[Ne];
        #else
        auto res = new double[Ne];
        #endif

        int runs = 1e7 / k + 1;
        double aggregate_time = 0;
        for (int p = 0; p < runs; ++p) {
            // Calculate
            double start = omp_get_wtime();
            spMV(ia, ja, a, b, Ne, res);
            double end = omp_get_wtime();

            aggregate_time += end - start;
        }

        double average_time = aggregate_time / runs;

        std::cout << 2 * ia[Ne] / (average_time * 1e9) << ",";
    }
    std::cout << std::endl;
    #endif

    return 0;
}
