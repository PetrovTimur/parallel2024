#include <omp.h>
#include <vector>
#include <thrust/device_vector.h>

#include "Solver/csr.h"
#include "Solver/solvers.cuh"
#include "Solver/Kernels/mathfunc.cuh"
#include "Utilities/input.h"

int main(int argc, char *argv[]) {
    omp_set_num_threads(omp_get_max_threads());

    int T  = omp_get_max_threads();
    std::cout << "T = " << T << std::endl;

    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

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

        auto a = new float[ia[Ne]];
        auto b = new float[Ne];
        auto diag = new float[Ne];

        fillCSR(ia, ja, a, b, diag, Ne);

        auto res = new float[Ne];

        int runs = 1e7 / k + 1;
        double aggregate_time = 0;
        int iterations = 0;
        for (int p = 0; p < runs; ++p) {
            // Calculate
            cudaEventRecord(start);
            iterations = solve(ia, ja, a, b, diag, Ne, res, 1e-3, 50);
            cudaEventRecord(stop);

            cudaEventSynchronize(stop);
            float milliseconds = 0;
            cudaEventElapsedTime(&milliseconds, start, stop);

            aggregate_time += milliseconds / 1000;
        }
        double average_time = aggregate_time / runs;

        std::cout << 2 * iterations * (0.5 * Ne + 2 * Ne + 3 * Ne + Ne + ia[Ne]) / (average_time * 1e9) << ",";
    }
    std::cout << std::endl;


}
