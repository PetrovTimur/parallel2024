#include <cmath>
#include <omp.h>
#include <iostream>
#include <ostream>
#include "Solver/Kernels/mathfunc.h"

#ifdef USE_MPI
#include <mpi.h>
#endif

int main(int argc, char *argv[]) {
    #ifdef USE_MPI
    omp_set_num_threads(omp_get_max_threads());

    std::vector<double> x;
    std::vector<double> y;
    double res, result;

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

        // FIll
        // #pragma omp parallel for proc_bind(spread)
        for (int i = 0; i < k / NumProc; i++) {
            x[i] = std::sin(MyID * k / NumProc + i);
            y[i] = std::cos(MyID * k / NumProc + i);
        }

        int runs = 1e9 / k + 1;
        double aggregate_time = 0;
        for (int p = 0; p < runs; ++p) {
            // Calculate
            double start = MPI_Wtime();
            dot(x.data(), y.data(), x.size(), res);
            MPI_Allreduce(&res, &result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            double end = MPI_Wtime();

            aggregate_time += end - start;

            // std::cout << end - start << std::endl;
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

    int T  = omp_get_max_threads();
    std::cout << "T = " << T << std::endl;

    for (int k = 1e6; k <= 1e8; k *= 10) {
        x.resize(k);
        y.resize(k);

        // FIll
        #pragma omp parallel default(none) shared(x, y, k)
        for (int i = 0; i < k; i++) {
            x[i] = std::sin(i);
            y[i] = std::cos(i);
        }

        int runs = 1e9 / k + 1;
        double aggregate_time = 0;
        for (int p = 0; p < runs; ++p) {
            double res;
            // Calculate
            double start = omp_get_wtime();
            dot(x.data(), y.data(), x.size(), res);
            double end = omp_get_wtime();

            aggregate_time += end - start;

            // std::cout << end - start << std::endl;
        }
        double average_time = aggregate_time / runs;

        std::cout << 2 * k / (average_time * 1e9) << ",";
    }
    std::cout << std::endl;
    #endif

    return 0;
}
