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

int main(int argc, char** argv) {
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
        std::vector<int> G2L((Nx + 1) * (Ny + 1), -1);
        std::vector<int> Part((Nx + 1) * (Ny + 1));

        std::tuple<int, int, int, int, int, int, int, int, int, int> t = input(Nx, Ny, Px, Py, MyID, L2G, G2L, Part);
        int i_start = std::get<0>(t);
        int i_end = std::get<1>(t);
        int i_count = std::get<2>(t);
        int j_start = std::get<3>(t);
        int j_end = std::get<4>(t);
        int j_count = std::get<5>(t);
        int top_halo = std::get<6>(t);
        int right_halo = std::get<7>(t);
        int bottom_halo = std::get<8>(t);
        int left_halo = std::get<9>(t);

        // std::vector<int> ia(N0 + 1);
        std::vector<int> ia = {0};
        // std::vector<int> ja(7 * (N0 + 1));
        std::vector<int> ja;

        makeCSR(Nx, Ny, K1, K2, i_start + top_halo, i_end - bottom_halo, j_start + left_halo, j_end - right_halo, G2L, ia, ja);

        std::vector<double> a(ja.size());
        std::vector<double> b(ia.size() - 1);
        std::vector<double> diag(ia.size() - 1);
        std::vector<double> res(ia.size() - 1);
        fillCSR(ia, ja, L2G, a, b, diag);

        std::vector<double> recv_buf;
        std::vector<int> recv_offset(7);

        std::vector<double> send_buf;
        std::vector<int> send_offset(7);

        ComInitOffsets(top_halo, left_halo, right_halo, bottom_halo, i_count, j_count, recv_offset, send_offset);
        recv_buf.resize(recv_offset[recv_offset.size() - 1]);
        send_buf.resize(send_offset[send_offset.size() - 1]);

        int runs = 1e7 / k + 1;
        double aggregate_time = 0;
        for (int p = 0; p < runs; ++p) {
            // Calculate
            double start = MPI_Wtime();
            Com(MyID, Px, top_halo, left_halo, right_halo, bottom_halo, i_count, j_count, b, recv_offset, send_offset, recv_buf, send_buf);
            spMV(ia, ja, a, b, recv_buf, res);
            MPI_Barrier(MPI_COMM_WORLD);
            double end = MPI_Wtime();

            aggregate_time += end - start;
        }

        double average_time = aggregate_time / runs;

        if (MyID == 0)
            std::cout << 2 * 6 * Nx * Ny / (average_time * 1e9) << ",";
    }
    if (MyID == 0)
        std::cout << std::endl;

    MPI_Finalize();

    #else
    omp_set_num_threads(omp_get_max_threads());
    std::vector<int> ia;
    std::vector<int> ja;
    std::vector<double> a;
    std::vector<double> b;
    std::vector<double> diag;
    std::vector<double> res;

    int T  = omp_get_max_threads();
    std::cout << "T = " << T << std::endl;

    for (int k = 1e5; k <= 1e8; k *= 10) {
        int K1 = 30;
        int K2 = 23;
        int Nx = 1000;
        int Ny = k / Nx / 5;

        auto stats = input(Nx, Ny, K1, K2);
        // std::tuple<int, int> t = input(Nx, Ny, K1, K2);
        // int nodes = std::get<0>(t);
        // int nonzero_elements = std::get<1>(t);
        int nodes = stats[1];
        int nonzero_elements = stats[4];

        ia.resize(nodes + 1);
        ja.resize(nonzero_elements);
        a.resize(nonzero_elements);
        b.resize(nodes);
        diag.resize(nodes);
        res.resize(nodes);

        makeCSR(Nx, Ny, K1, K2, ia, ja);
        fillCSR(ia.data(), ja.data(), a.data(), b.data(), diag.data(), nodes);

        int runs = 1e7 / k + 1;
        double aggregate_time = 0;
        for (int p = 0; p < runs; ++p) {
            // Calculate
            double start = omp_get_wtime();
            spMV(ia, ja, a, b, res);
            double end = omp_get_wtime();

            aggregate_time += end - start;
        }

        double average_time = aggregate_time / runs;

        std::cout << 2 * a.size() / (average_time * 1e9) << ", ";

        delete[] stats;
    }
    std::cout << std::endl;
    #endif

    return 0;
}
