#include <cassert>
#include <cmath>
#include <omp.h>
#include <iostream>
#include <ostream>
#include <vector>

#include "Solver/Kernels/mathfunc.h"

#ifdef USE_MPI
#include <mpi.h>
#endif

int main() {
    #ifdef USE_MPI
    omp_set_num_threads(omp_get_max_threads());

    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> res;

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
        x.resize(k / NumProc);
        y.resize(k / NumProc);
        res.resize(k / NumProc);

        // Fill
        // #pragma omp parallel for proc_bind(spread)
        for (int i = 0; i < k / NumProc; i++) {
            x[i] = sin(MyID * k / NumProc + i);
            y[i] = cos(MyID * k / NumProc + i);
        }

        int runs = 1e9 / k + 1;
        double aggregate_time = 0;
        for (int p = 1; p < runs + 1; ++p) {
            double alpha = std::tan(p);
            // Calculate
            double start = MPI_Wtime();
            axpy(alpha, x, y, res);
            MPI_Barrier(MPI_COMM_WORLD);
            double end = MPI_Wtime();

            aggregate_time += end - start;

        }

        double average_time = aggregate_time / runs;

        if (MyID == 0)
            std::cout << 2 * k / (average_time * 1e9) << ",";
    }
    if (MyID == 0)
        std::cout << std::endl;

    MPI_Finalize();

    #else
    omp_set_num_threads(omp_get_max_threads());

    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> res;

    int T  = omp_get_max_threads();
    std::cout << "T = " << T << std::endl;

    for (int k = 1e5; k <= 1e8; k *= 10) {
        x.resize(k);
        y.resize(k);
        res.resize(k);

        // Fill
        #pragma omp parallel for proc_bind(spread)
        for (int i = 0; i < k; i++) {
            x[i] = sin(i);
            y[i] = cos(i);
        }

        int runs = 1e9 / k + 1;
        double aggregate_time = 0;
        for (int p = 1; p < runs + 1; ++p) {
            double alpha = std::tan(p);
            // Calculate
            double start = omp_get_wtime();
            axpy(alpha, x.data(), y.data(), x.size(), res.data());
            double end = omp_get_wtime();

            aggregate_time += end - start;

        }

        double average_time = aggregate_time / runs;

        std::cout << 2 * k / (average_time * 1e9) << ", ";
        std::cout << average_time << std::endl;
    }
    std::cout << std::endl;
    #endif

    return 0;
}
