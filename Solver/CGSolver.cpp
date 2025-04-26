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



int main(int argc, char** argv) {
    arguments args{};

    args.log_dir = "";
    args.eps = 1e-3;
    args.maxit = 100;

    argp_parse (&argp, argc, argv, 0, nullptr, &args);

    if (!args.log_dir.empty()) {
        Logger::setLogDirectory(args.log_dir);
        LOG_INFO << "Starting with log directory: " << args.log_dir << std::endl;
    }

    omp_set_num_threads(omp_get_max_threads());

    int Nx = args.Nx;
    checkInput(Nx, "Nx");
    int Ny = args.Ny;
    checkInput(Ny, "Ny");
    int K1 = args.K1;
    checkInput(K1, "K1");
    int K2 = args.K2;
    checkInput(K2, "K2");

    double eps = args.eps;
    checkInput(eps, "eps");
    int maxit = args.maxit;
    checkInput(maxit, "maxit");

    LOG_INFO << "Nx = " << Nx << ", Ny = " << Ny << ", K1 = " << K1 << ", K2 = " << K2 << std::endl;

    const gridInfo grid(Nx, Ny, K1, K2);

    int Nn = grid.totalNodes;
    int Ne = grid.totalElements;
    // int nnz = 4 * stats[0] + 2 * stats[2];

    LOG_INFO << "Nn = " << Nn << ", Ne = " << Ne << std::endl;

    LOG_INFO << "Using " << omp_get_max_threads() << " threads" << std::endl;

    std::vector<int> ia_en;
    std::vector<int> ja_en;

    double start = omp_get_wtime();
    makeIncidenceMatrixCSR(Nx, Ny, K1, K2, 0, Ny - 1, 0, Nx - 1, ia_en, ja_en);
    double end = omp_get_wtime();

    LOG_INFO << "EN construction done in " << end - start << " seconds" << std::endl;

    // int *ia_en = std::get<0>(matrix);
    // int *ja_en = std::get<1>(matrix);

    #ifdef DEBUG_MODE
    LOG_DEBUG << "IA_EN:\t";
    printVector(ia_en, LOG);
    LOG_DEBUG << "JA_EN:\t";
    printVector(ja_en, LOG);
    #endif

    // std::vector<int> ia_ne;
    // std::vector<int> ja_ne;

    int *ia_ne, *ja_ne;

    start = omp_get_wtime();
    transposeCSR(ia_en, ja_en, Nn, ia_ne, ja_ne);
    end = omp_get_wtime();

    LOG_INFO << "Transpose (NE) construction done in " << end - start << " seconds" << std::endl;

    // int *ia_ne = std::get<0>(matrix_transposed);
    // int *ja_ne = std::get<1>(matrix_transposed);

    #ifdef DEBUG_MODE
    LOG_DEBUG << "IA_NE:\t";
    printArray(ia_ne, Nn + 1, LOG);
    LOG_DEBUG << "JA_NE:\t";
    printArray(ja_ne, ia_ne[Nn], LOG);
    #endif

    start = omp_get_wtime();
    auto matrix_ee = buildAdjacencyMatrixCSRUsingSort(ia_en.data(), ja_en.data(), ia_ne, ja_ne, Ne, Nn);
    end = omp_get_wtime();

    LOG_INFO << "Adjacency matrix (EE) construction done in " << end - start << " seconds" << std::endl;

    // Ne = Nn;
    // auto matrix_nn = buildAdjacencyMatrixCSRUsingSort(ia_ne.data(), ja_ne.data(), ia_en.data(), ja_en.data(), Ne, Nn);
    int *ia_ee = std::get<0>(matrix_ee);
    int *ja_ee = std::get<1>(matrix_ee);

    // std::cout << "Built adj matrix" << std::endl;

    int *ia = ia_ee;
    int *ja = ja_ee;

    #ifdef DEBUG_MODE
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

    start = omp_get_wtime();
    fillCSR(ia, ja, a, b, diag, Ne);
    end = omp_get_wtime();

    LOG_INFO << "Fill done in " << end - start << " seconds" << std::endl;

    #ifdef DEBUG_MODE
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

    start = omp_get_wtime();
    int iterations = solve(ia, ja, a, b, diag, Ne, res, eps, maxit);
    end = omp_get_wtime();

    LOG_INFO << "Work done in " << end - start << " seconds" << std::endl;
    LOG_INFO << "Convergence required "  << iterations << " iterations" << std::endl;
    LOG_INFO << res[0] << std::endl;

    // delete[] ia_en;
    // delete[] ja_en;
    delete[] ia_ne;
    delete[] ja_ne;
    delete[] ia;
    delete[] ja;
    delete[] a;
    delete[] b;
    delete[] diag;

    #ifdef DEBUG_MODE
    LOG_DEBUG << "res:\t";
    printArray(res, Ne, LOG);
    #endif

    delete[] res;

    return 0;
}
