#include <unordered_map>

#include "omp.h"
#include "csr.h"
#include "Utilities/input.h"
#include "Utilities/logger.h"
#include "Utilities/argparse.h"

#ifdef USE_CUDA
#include "solvers.cuh"
#else
#include "solvers.h"
#endif

#ifdef USE_MPI
#include <mpi.h>
#include "Utilities/coms.h"
#endif


int main(int argc, char** argv) {
    #ifdef USE_MPI
    struct arguments arguments{};

    /* Default values. */
    arguments.output_file = "-";
    arguments.eps = 1e-3;
    arguments.maxit = 1000;

    /* Parse our arguments; every option seen by parse_opt will
       be reflected in arguments. */
    argp_parse (&argp, argc, argv, 0, nullptr, &arguments);

    omp_set_num_threads(1);

    const int Nx = arguments.Nx;
    const int Ny = arguments.Ny;
    const int K1 = arguments.K1;
    const int K2 = arguments.K2;
    const int Px = arguments.Px;
    const int Py = arguments.Py;

    double eps = arguments.eps;
    int maxit = arguments.maxit;

    int mpi_res = MPI_Init(&argc, &argv); // первым делом подключаемся к MPI
    if(mpi_res != MPI_SUCCESS) {
        MPI_Abort(MPI_COMM_WORLD, -1);
        exit(mpi_res);
    }

    int NumProc;
    mpi_res = MPI_Comm_size(MPI_COMM_WORLD,&NumProc); // узнаем число процессов
    if(mpi_res != MPI_SUCCESS){
        MPI_Abort(MPI_COMM_WORLD, -1);
        exit(mpi_res);
    }

    int MyID;
    mpi_res = MPI_Comm_rank(MPI_COMM_WORLD,&MyID); // узнаем номер данного процесса
    if(mpi_res != MPI_SUCCESS){
        MPI_Abort(MPI_COMM_WORLD, -1);
        exit(mpi_res);
    }


    if (MyID == 0) {
        LOG_INFO << "Nx = " << Nx << ", Ny = " << Ny << ", K1 = " << K1 << ", K2 = " << K2 << ", Px = " << Px << ", Py = " << Py << std::endl;
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
    // std::vector<int> G2L;
    std::vector<int> Part;

    auto t = input(Nx, Ny, Px, Py, MyID, L2G, Part);
    // int i_start = std::get<0>(t);
    // int i_end = std::get<1>(t);
    // int i_count = std::get<2>(t);
    // int j_start = std::get<3>(t);
    // int j_end = std::get<4>(t);
    // int j_count = std::get<5>(t);
    // int top_halo = std::get<6>(t);
    // int right_halo = std::get<7>(t);
    // int bottom_halo = std::get<8>(t);
    // int left_halo = std::get<9>(t);
    // int N_local = std::get<10>(t);


    // std::cout << "My ID = " << MyID << ", i_start " << i_start << ", i_count " << i_count << ", i_end " << i_end << ", j_start " << j_start << ", j_count " << j_count << ", j_end " << j_end << std::endl;
    // std::cout << "My ID = " << MyID << ", top_halo " << top_halo << ", right_halo " << right_halo << ", bottom_halo " << bottom_halo << ", left_halo " << left_halo << std::endl;


    // std::vector<int> ia(N0 + 1);
    std::vector<int> ia_en;
    // std::vector<int> ja(7 * (N0 + 1));
    std::vector<int> ja_en;

    // std::cout << "My ID = " << MyID << std::endl;
    // std::cout << "L2G size = " << L2G.size() << std::endl;
    // printVector(L2G, std::cout);
    //
    // MPI_Barrier(MPI_COMM_WORLD);

    makeIncidenceMatrixCSR(Nx, Ny, K1, K2, L2G, ia_en, ja_en, Part);


    std::unordered_map<int, int> G2L;
    fillG2L(L2G, G2L);

    std::unordered_map<int, int> G2L_nodes;
    constructG2L(ia_en, ja_en, G2L_nodes);

    localizeCSR(ia_en.data(), ia_en.size(), ja_en.data(), G2L_nodes);

    // printVector(L2G, std::cout);
    // std::cout << "G2L_nodes size = " << G2L_nodes.size() << std::endl;
    // printVector(ja_en, std::cout);

    std::vector<int> ia_ne, ja_ne;
    transposeCSR(ia_en, ja_en, G2L_nodes.size(), ia_ne, ja_ne);

    // printVector(ja_ne, std::cout);

    // MPI_Barrier(MPI_COMM_WORLD);
    // std::cout << "My ID = " << MyID << ", Part size " << Part.size() << ", Ne = " << ia_en.size() - 1 << std::endl;
    // printVector(Part, std::cout);

    std::vector<int> ia_ee, ja_ee;
    buildAdjacencyMatrixCSRUsingSort(ia_en.data(), ja_en.data(), ia_ne.data(), ja_ne.data(), ia_ee, ja_ee, ia_en.size() - 1, Part, MyID);
    // auto ia_ee = std::get<0>(matrix);
    // auto ja_ee = std::get<1>(matrix);

    // MPI_Barrier(MPI_COMM_WORLD);
    // std::cout << "My ID = " << MyID << std::endl;
    // printVector(L2G, std::cout);

    std::vector<double> a(ja_ee.size());
    std::vector<double> b(ia_ee.size() - 1);
    std::vector<double> diag(ia_ee.size() - 1);
    fillCSR(ia_ee, ja_ee, L2G, a, b, diag);

    // std::cout << "My ID = " << MyID << std::endl;
    // printVector(ja_ee, std::cout);

    // printVector(ia);
    // printVector(ja);
    // if (MyID == 0) {
    //     LOG_INFO << "My ID = " << MyID << ", b size: " << b.size() << std::endl;
    //     printVector(ia_ee, LOG);
    //     printVector(a, LOG);
    //     printVector(b, LOG);
    // }
    // printVector(b);

    std::vector<double> res(ia_ee.size() - 1);

    double start = MPI_Wtime();
    int iterations = solve(MyID, Part, L2G, ia_ee, ja_ee, a, b, diag, res);


    // std::cout << "MyID: " << MyID << ", x[0]: " << res[0] << std::endl;
    // for (double x : res)
    //     std::cout << x << " ";
    // std::cout << std::endl;

    // out.close();
    MPI_Finalize();
    double end = MPI_Wtime();

    if (MyID == 0) {
        LOG_INFO << "Work took " << end - start << " seconds" << std::endl;
        LOG_INFO << "Convergence required "  << iterations << " iterations" << std::endl;
        LOG_INFO << res[0] << std::endl;
    }

    #else
    struct arguments arguments{};

    /* Default values. */
    arguments.output_file = "-";
    arguments.eps = 1e-3;
    arguments.maxit = 100;

    argp_parse (&argp, argc, argv, 0, nullptr, &arguments);

    omp_set_num_threads(omp_get_max_threads());

    int Nx = arguments.Nx;
    int Ny = arguments.Ny;
    int K1 = arguments.K1;
    int K2 = arguments.K2;

    double eps = arguments.eps;
    int maxit = arguments.maxit;

    LOG_INFO << "Nx = " << Nx << ", Ny = " << Ny << ", K1 = " << K1 << ", K2 = " << K2 << std::endl;

    auto stats = input(Nx, Ny, K1, K2);
    int nodes = stats[1];

    int Nn = nodes;
    int Ne = stats[0] + stats[2];
    int nnz = 4 * stats[0] + 2 * stats[2];

    LOG_DEBUG << "Nn = " << Nn << ", Ne = " << Ne << ", nnz = " << nnz << std::endl;

    std::vector<int> ia_en;
    std::vector<int> ja_en;
    makeIncidenceMatrixCSR(Nx, Ny, K1, K2, 0, Ny - 1, 0, Nx - 1, ia_en, ja_en);
    // int *ia_en = std::get<0>(matrix);
    // int *ja_en = std::get<1>(matrix);

    #ifdef USE_DEBUG_MODE
    LOG_DEBUG << "IA_EN:\t";
    printVector(ia_en, LOG);
    LOG_DEBUG << "JA_EN:\t";
    printVector(ja_en, LOG);
    #endif

    std::vector<int> ia_ne;
    std::vector<int> ja_ne;
    transposeCSR(ia_en, ja_en, Nn, ia_ne, ja_ne);
    // int *ia_ne = std::get<0>(matrix_transposed);
    // int *ja_ne = std::get<1>(matrix_transposed);

    #ifdef USE_DEBUG_MODE
    LOG_DEBUG << "IA_NE:\t";
    printArray(ia_ne.data(), Nn + 1, LOG);
    LOG_DEBUG << "JA_NE:\t";
    printArray(ja_ne.data(), ia_ne[Nn], LOG);
    #endif

    auto matrix_nn = buildAdjacencyMatrixCSRUsingSort(ia_en.data(), ja_en.data(), ia_ne.data(), ja_ne.data(), Ne, Nn);
    // Ne = Nn;
    // auto matrix_nn = buildAdjacencyMatrixCSRUsingSort(ia_ne.data(), ja_ne.data(), ia_en.data(), ja_en.data(), Ne, Nn);
    int *ia_nn = std::get<0>(matrix_nn);
    int *ja_nn = std::get<1>(matrix_nn);

    std::cout << "Built adj matrix" << std::endl;

    int *ia = ia_nn;
    int *ja = ja_nn;

    #ifdef USE_DEBUG_MODE
    // std::vector<std::vector<bool>> matrix(nodes + 1, std::vector<bool>(nodes + 1, false));
    // buildMatrixFromCSR(ia, ja, matrix);
    // LOG_DEBUG << "Matrix:" << std::endl;
    // printMatrix(matrix, LOG);
    #endif

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

    #ifdef USE_DEBUG_MODE
    LOG_DEBUG << "IA:\t\t";
    printArray(ia, Ne + 1, LOG);

    LOG_DEBUG << "JA:\t\t";
    printArray(ja, ia[Ne],LOG);

    LOG_DEBUG << "A:\t\t";
    printArray(a, ia[Ne], LOG);

    LOG_DEBUG << "b:\t\t";
    printArray(b, Ne, LOG);
    #endif

    LOG << std::endl;

    #ifdef USE_CUDA
    auto res = new float[Ne];
    #else
    auto res = new double[Ne];
    #endif

    double start = omp_get_wtime();
    int iterations = solve(ia, ja, a, b, diag, Ne + 1, res, eps, maxit);
    double end = omp_get_wtime();

    LOG_INFO << "Work took " << end - start << " seconds" << std::endl;
    LOG_INFO << "Convergence required "  << iterations << " iterations" << std::endl;
    LOG_INFO << res[0] << std::endl;

    delete[] stats;
    // delete[] ia_en;
    // delete[] ja_en;
    // delete[] ia_ne;
    // delete[] ja_ne;
    delete[] ia;
    delete[] ja;
    delete[] a;
    delete[] b;
    delete[] diag;

    #ifdef USE_DEBUG_MODE
    LOG_DEBUG << "res:\t";
    printArray(res, Ne, LOG);
    #endif

    delete[] res;

    #endif

    return 0;
}
