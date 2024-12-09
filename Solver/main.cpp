#include <iostream>
#include <string>
#include <tuple>
#include <vector>
#include <fstream>
#include <omp.h>
#include <unistd.h>

#include "csr.h"
#include "Utilities/input.h"
#include "solvers.h"
#include "Kernels/mathfunc.h"
#include "Utilities/argparse.h"
#ifdef USE_MPI
#include "Utilities/coms.h"
#endif

#ifdef USE_MPI
#include <mpi.h>
#endif


int main(int argc, char** argv) {
    #ifdef USE_MPI
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

    if (MyID == 0) {
        std::cout << "Nx = " << Nx << std::endl;
        std::cout << "Ny = " << Ny << std::endl;
        std::cout << "K1 = " << K1 << std::endl;
        std::cout << "K2 = " << K2 << std::endl;
        std::cout << "Px = " << Px << std::endl;
        std::cout << "Py = " << Py << std::endl;
    }

    int MyID_j = MyID % Px;
    int MyID_i = MyID / Px;

    int i_start, i_end, i_count, j_start, j_end, j_count;
    j_start = MyID_j * ((Nx + 1) / Px) + std::min((Nx + 1) % Px, MyID_j);
    j_count = (Nx + 1) / Px + ((Nx + 1) % Px > MyID_j);
    j_end = j_start + j_count - 1;
    int left_halo = MyID_j > 0; // Halo
    j_start -= left_halo;
    int right_halo = MyID_j < Px - 1; // Halo
    j_end += right_halo;

    i_start = MyID_i * ((Ny + 1) / Py) + std::min((Ny + 1) % Py, MyID_i);
    i_count = (Ny + 1) / Py + ((Ny + 1) % Py > MyID_i);
    i_end = i_start + i_count - 1;
    int top_halo = MyID_i > 0; // Halo
    i_start -= top_halo;
    int bottom_halo = MyID_i < Py - 1; // Halo
    i_end += bottom_halo;

    // std::cout << "MyID: " << MyID << std::endl;
    // std::cout << "i_start: " << i_start << ", i_end: " << i_end << std::endl;
    // std::cout << "j_start: " << j_start << ", j_end: " << j_end << std::endl;

    int k = 0, N0;
    std::vector<int> L2G((i_end - i_start + 1) * (j_end - j_start + 1));
    for (int i = i_start + top_halo; i <= i_end - bottom_halo; i++) {
        for (int j = j_start + left_halo; j <= j_end - right_halo; j++) {
            L2G[k++] = i * (Nx + 1) + j;
        }
    }
    N0 = k;
    if (top_halo) {
        for (int j = j_start + left_halo; j <= j_end - right_halo; j++) {
            L2G[k++] = i_start * (Nx + 1) + j;
        }
        if (right_halo)
            L2G[k++] = i_start * (Nx + 1) + j_end;
    }
    if (left_halo) {
        for (int i = i_start + top_halo; i <= i_end - bottom_halo; i++) {
            L2G[k++] = i * (Nx + 1) + j_start;
        }
    }
    if (right_halo) {
        for (int i = i_start + top_halo; i <= i_end - bottom_halo; i++) {
            L2G[k++] = i * (Nx + 1) + j_end;
        }
    }
    if (bottom_halo) {
        if (left_halo)
            L2G[k++] = i_end * (Nx + 1) + j_start;
        for (int j = j_start + left_halo; j <= j_end - right_halo; j++) {
            L2G[k++] = i_end * (Nx + 1) + j;
        }
    }

    L2G.resize(k);

    // if (MyID == 7) {
    //     std::cout << i_start << i_end << j_start << j_end << std::endl;
    //     // std::cout << top_halo << bottom_halo << left_halo << right_halo << std::endl;
    //     std::cout << "L2G size: " << L2G.size() << std::endl;
    //     for (int i = 0; i < L2G.size(); i++) {
    //         std::cout << L2G[i] << " ";
    //     }
    //     std::cout << std::endl;
    // }

    std::vector<int> G2L((Nx + 1) * (Ny + 1), -1);
    for (int iL = 0; iL < L2G.size(); iL++) {
        G2L[L2G[iL]] = iL;
    }


    std::vector<int> Part((Nx + 1) * (Ny + 1));
    for (int i = 0; i < Py; i++) {
        int is = i * ((Ny + 1) / Py) + std::min((Ny + 1) % Py, i);
        int ic = (Ny + 1) / Py + ((Ny + 1) % Py > i) - 1;
        int ie = is + ic;
        for (int j = 0; j < Px; j++) {
            int js = j * ((Nx + 1) / Px) + std::min((Nx + 1) % Px, j);
            int jc = (Nx + 1) / Px + ((Nx + 1) % Px > j) - 1;
            int je = js + jc;

#ifdef  USE_DEBUG_MODE
            if (MyID == 0) {
                std::cout << "id: " << i * Px + j << std::endl;
                std::cout << is << " " << ie << " " << js << " " << jc << std::endl;
            }
#endif

            for (int ii = is; ii <= ie; ii++) {
                for (int jj = js; jj <= je; jj++) {
                    Part[ii * (Nx + 1) + jj] = i * Px + j;
                }
            }
        }
    }

#ifdef USE_DEBUG_MODE
    if (MyID == 0) {
        std::cout << "Part size: " << Part.size() << std::endl;
        for (int i : Part) {
            std::cout << i << " ";
        }
        std::cout << std::endl;
    }
#endif

    std::vector<int> ia(N0 + 1);
    std::vector<int> ja(7 * (N0 + 1));

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
    solve(MyID, Px, top_halo, left_halo, right_halo, bottom_halo, i_count, j_count, ia, ja, a, b, diag, res);

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

    std::fstream out;
    if (arguments.output_file[0] != '-') {
        out.open(arguments.output_file, std::ios::out);
        std::cout.rdbuf(out.rdbuf()); //save and redirect
    }

    int Nx = std::stoi(arguments.args[0]);
    int Ny = std::stoi(arguments.args[1]);
    int K1 = std::stoi(arguments.args[2]);
    int K2 = std::stoi(arguments.args[3]);

    std::cout << "Nx = " << Nx << std::endl;
    std::cout << "Ny = " << Ny << std::endl;
    std::cout << "K1 = " << K1 << std::endl;
    std::cout << "K2 = " << K2 << std::endl;

    std::tuple<int, int> t = input(Nx, Ny, K1, K2);
    int nodes = std::get<0>(t);
    int nonzero_elements =std::get<1>(t);

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

    fillCSR(ia, ja, a, b, diag);

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

    int iterations = solve(ia, ja, a, b, diag, res);

    printVector(res);

    double end = omp_get_wtime();
    std::cout << "Work took " << end - start << " seconds\n";
    std::cout << "Convergence required "  << iterations << " iterations\n";

    out.close();

    #ifdef USE_DEBUG_MODE
    std::cout << "res: ";
    printVector(res);
    #endif

    #endif

    return 0;
}
