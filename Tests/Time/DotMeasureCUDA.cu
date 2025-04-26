#include <omp.h>
#include <vector>
#include <thrust/device_vector.h>

#include "Solver/Kernels/mathfunc.cuh"

int main(int argc, char *argv[]) {
    omp_set_num_threads(omp_get_max_threads());

    int threads = 256;
    int blocks = 204;

    std::vector<float> x;
    std::vector<float> y;

    int T  = omp_get_max_threads();
    std::cout << "T = " << T << std::endl;

    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    for (int k = 1e6; k <= 1e8; k *= 10) {
        x.resize(k);
        y.resize(k);

        // FIll
        #pragma omp parallel default(none) shared(x, y, k)
        for (int i = 0; i < k; i++) {
            x[i] = std::sin(i);
            y[i] = std::cos(i);
        }

        auto d_x = thrust::device_vector<float>(x.begin(), x.end());
        auto d_y = thrust::device_vector<float>(y.begin(), y.end());
        auto buf = thrust::device_vector<float> (blocks);

        float res_host;
        float *res;
        cudaMalloc(&res, sizeof(float));

        int runs = 1e9 / k + 1;
        double aggregate_time = 0;
        for (int p = 0; p < runs; ++p) {
            // Calculate
            cudaEventRecord(start);
            dot_gpu(threads, blocks, d_x.data().get(), d_y.data().get(), buf.data().get(), res, k);
            cudaMemcpy(&res_host, res, sizeof(float), cudaMemcpyDeviceToHost);
            cudaEventRecord(stop);

            cudaEventSynchronize(stop);
            float milliseconds = 0;
            cudaEventElapsedTime(&milliseconds, start, stop);
            // std::cout << 2 * k / 1e9 / (milliseconds / 1000.0) << " GFLOPS, " << milliseconds << " ms" << std::endl;

            aggregate_time += milliseconds / 1000;

            // std::cout << end - start << std::endl;
        }
        double average_time = aggregate_time / runs;

        std::cout << 2 * k / (average_time * 1e9) << ",";
    }
    std::cout << std::endl;


}
