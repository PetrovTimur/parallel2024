#include <iostream>
#include <string>
#include <tuple>
#include <vector>
#include <fstream>
#include <omp.h>
#include <unistd.h>

#include "csr.h"
#include "Utilities/input.h"
#include "Utilities/logger.h"
#include "solvers.h"
#include "Kernels/mathfunc.h"
#include "Utilities/argparse.h"

#ifdef USE_MPI
#include <mpi.h>
#include "Utilities/coms.h"
#endif


int main(int argc, char** argv) {
    #ifdef USE_MPI
    struct arguments arguments{};

    /* Default values. */
    arguments.output_file = "-";

    /* Parse our arguments; every option seen by parse_opt will
       be reflected in arguments. */
    argp_parse (&argp, argc, argv, 0, nullptr, &arguments);

    // omp_set_num_threads(omp_get_max_threads());

    std::fstream out;
    if (arguments.output_file[0] != '-') {
        out.open(arguments.output_file, std::ios::out);
        std::cout.rdbuf(out.rdbuf()); //save and redirect
    }

    int Nx = std::stoi(arguments.args[0]);
    int Ny = std::stoi(arguments.args[1]);
    int K1 = std::stoi(arguments.args[2]);
    int K2 = std::stoi(arguments.args[3]);
    int Px = std::stoi(arguments.args[4]);
    int Py = std::stoi(arguments.args[5]);

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
        std::cout << "Nx = " << Nx << std::endl;
        std::cout << "Ny = " << Ny << std::endl;
        std::cout << "K1 = " << K1 << std::endl;
        std::cout << "K2 = " << K2 << std::endl;
        std::cout << "Px = " << Px << std::endl;
        std::cout << "Py = " << Py << std::endl;
        std::cout << "M = " << omp_get_max_threads() << std::endl;
    }

    if (NumProc != Px * Py) {
        if (MyID == 0) {
            std::cout << "Number of processes Px*Py doesn't match MPI config" << std::endl;
        }

        MPI_Finalize();
        return 0;
    }

    // std::cout << omp_get_max_threads() << std::endl;
    // omp_set_num_threads(omp_get_max_threads() / P0);
    // std::cout << omp_get_max_threads() << std::endl;

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
    fillCSR(ia, ja, L2G, a, b, diag);

    // printVector(ia);
    // printVector(ja);
    // printVector(a);
    // printVector(b);

    std::vector<double> res(b.size());

    double start = MPI_Wtime();
    int iterations = solve(MyID, Px, top_halo, left_halo, right_halo, bottom_halo, i_count, j_count, ia, ja, a, b, diag, res);
    std::cout << "Iterations: " << iterations << std::endl;

    std::cout << "MyID: " << MyID << ", x[0]: " << res[0] << std::endl;
    // for (double x : res)
    //     std::cout << x << " ";
    // std::cout << std::endl;

    out.close();
    MPI_Finalize();

    if (MyID == 0) {
        double end = MPI_Wtime();
        std::cout << end - start << " sec" << std::endl;
    }

    #else
    struct arguments arguments{};

    /* Default values. */
    arguments.output_file = "-";

    /* Parse our arguments; every option seen by parse_opt will
       be reflected in arguments. */
    argp_parse (&argp, argc, argv, 0, nullptr, &arguments);

    omp_set_num_threads(omp_get_max_threads());

    int Nx = std::stoi(arguments.args[0]);
    int Ny = std::stoi(arguments.args[1]);
    int K1 = std::stoi(arguments.args[2]);
    int K2 = std::stoi(arguments.args[3]);

    LOG_INFO << "Nx = " << Nx << ", Ny = " << Ny << ", K1 = " << K1 << ", K2 = " << K2 << std::endl;

    auto stats = input(Nx, Ny, K1, K2);
    int nodes = stats[1];
    int nonzero_elements = stats[4];

    std::vector<int> ia(nodes + 1);
    std::vector<int> ja(nonzero_elements);
    std::vector<double> a(nonzero_elements);
    std::vector<double> b(nodes);
    std::vector<double> diag(nodes);

    makeCSR(Nx, Ny, K1, K2, ia, ja);

    #ifdef USE_DEBUG_MODE
    std::vector<std::vector<bool>> matrix(nodes + 1, std::vector<bool>(nodes + 1, false));
    buildAdjacencyMatrix(ia, ja, matrix);
    printMatrix(matrix);
    #endif

    fillCSR(ia.data(), ja.data(), a.data(), b.data(), diag.data(), nodes);

    #ifdef USE_DEBUG_MODE
    std::cout << "IA: ";
    printVector(ia);

    std::cout << "JA: ";
    printVector(ja);

    std::cout << "A: ";
    printVector(a);

    std::cout << "b: ";
    printVector(b);
    #endif
    std::vector<double> res(nodes);

    double start = omp_get_wtime();
    int iterations = solve(ia.data(), ja.data(), a.data(), b.data(), diag.data(), ia.size(), res.data());
    double end = omp_get_wtime();

    printVector(res);
    LOG_INFO << "Work took " << end - start << " seconds" << std::endl;
    LOG_INFO << "Convergence required "  << iterations << " iterations" << std::endl;

    // out.close();
    delete[] stats;

    #ifdef USE_DEBUG_MODE
    std::cout << "res: ";
    printVector(res);
    #endif

    #endif

    return 0;
}
