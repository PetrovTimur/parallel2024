#include <omp.h>
#include <vector>
#include <thrust/device_vector.h>

#include "Solver/csr.h"
#include "Solver/Kernels/mathfunc.cuh"
#include "Utilities/input.h"

int main(int argc, char *argv[]) {
    omp_set_num_threads(omp_get_max_threads());

    int threads = 256;
    int blocks = 204;

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

        thrust::device_vector<int> d_ia(ia, ia + Ne + 1);
        thrust::device_vector<int> d_ja(ja, ja + ia[Ne]);
        thrust::device_vector<float> d_a(a, a + ia[Ne]);
        thrust::device_vector<float> d_b(b, b + Ne);
        thrust::device_vector<float> d_diag(diag, diag + Ne);
        thrust::device_vector<float> res(Ne);


        int runs = 1e7 / k + 1;
        double aggregate_time = 0;
        for (int p = 0; p < runs; ++p) {
            // Calculate
            cudaEventRecord(start);
            spMV<<<blocks, threads>>>(d_ia.data().get(), d_ja.data().get(), d_a.data().get(), d_b.data().get(), res.data().get(), Ne);
            cudaEventRecord(stop);

            cudaEventSynchronize(stop);
            float milliseconds = 0;
            cudaEventElapsedTime(&milliseconds, start, stop);

            aggregate_time += milliseconds / 1000;
        }
        double average_time = aggregate_time / runs;

        std::cout << 2 * k / (average_time * 1e9) << ",";
    }
    std::cout << std::endl;


}
