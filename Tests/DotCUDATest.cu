#include "thrust/device_vector.h"
#include <thrust/host_vector.h>
#include "Solver/Kernels/mathfunc.cuh"
#include "Utilities/cuda_helper.cuh"
#include "cuda_runtime.h"
#include <iostream>


int main(int argc, char *argv[]) {
    int blocks, threads;
    getDeviceSpecs(blocks, threads);

    if (argc > 2) {
        blocks = std::atoi(argv[1]);
        threads = std::atoi(argv[2]);
    }

    std::cout << blocks << " blocks, " << threads << " threads" << std::endl;

    int size = 100000000;

    if (argc > 3) {
        size = std::atoi(argv[3]);
    }

    std::cout << size << " numbers" << std::endl;

    float a = 2.f;
    float b = 3.f;

    thrust::device_vector<float> d_x(size, a);
    thrust::device_vector<float> d_y(size, b);
    thrust::device_vector<float> d_z(blocks);

    float res_host;
    float *res;
    cudaMalloc(&res, sizeof(float));

    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    cudaEventRecord(start);
    dot_gpu(threads, blocks, d_x.data().get(), d_y.data().get(), d_z.data().get(), res, size);
    cudaEventRecord(stop);

    cudaEventSynchronize(stop);
    float milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, start, stop);
    std::cout << 2 * size / 1e9 / (milliseconds / 1000.0) << " GFLOPS, " << milliseconds / 1000.0 << " s" << std::endl;


    cudaEventDestroy(start);
    cudaEventDestroy(stop);

    cudaMemcpy(&res_host, res, sizeof(float), cudaMemcpyDeviceToHost);

    // Verify the result & print first error if any
    bool correct = (res_host == a * b * size);

    std::cout << res_host << std::endl;
    std::cout << size * a * b << std::endl;


    if (correct) {
        std::cout << "Dot test passed!" << std::endl;
    } else {
        std::cout << "Dot test failed!" << std::endl;
    }

    cudaFree(res);

    return 0;
}