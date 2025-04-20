#include "thrust/device_vector.h"
#include <thrust/host_vector.h>
#include "Solver/Kernels/mathfunc.cuh"
#include "cuda_runtime.h"
#include <iostream>
#include <omp.h>
#include <vector>

int main() {
    int blocks = 128;
    int threads = 256;
    int size = 100000000;
    float a = 1.0f;

    thrust::host_vector<float> x(size, a);
    thrust::device_vector<float> d_x(size);
    thrust::device_vector<float> d_y(blocks);
    d_x = x;

    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    cudaEventRecord(start);
    reduce0<<<blocks, threads, threads*sizeof(float)>>>(d_x.data().get(), d_y.data().get(), size);
    reduce0<<<1, blocks, blocks*sizeof(float)>>>(d_y.data().get(), d_x.data().get(), blocks);
    cudaEventRecord(stop);

    cudaEventSynchronize(stop);
    float milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, start, stop);
    std::cout << 2 * size / 1e9 / (milliseconds / 1000.0) << " GFLOPS, " << milliseconds << " ms" << std::endl;

    float res = d_x[0];
    d_x[0] = a;

    // Verify the result & print first error if any
    bool correct = (res == size * a);

    std::cout << res << std::endl;
    std::cout << size * a << std::endl;


    if (correct) {
        std::cout << "Reduce test passed!" << std::endl;
    } else {
        std::cout << "Reduce test failed!" << std::endl;
    }

    std::cout << std::endl << "-----------------" << std::endl << std::endl;


    blocks = (size + threads - 1) / threads;
    thrust::device_vector<float> d_intermediate(blocks);
    d_y.resize(blocks);

    cudaEventRecord(start);
    reduce1<<<blocks, threads>>>(d_x.data().get(), d_y.data().get(), size);

    // std::cout << d_y[0] << d_y[1] << std::endl;

    int s = blocks;
    while (s > 1) {
        blocks = (s + threads - 1) / threads;

        cudaMemcpy(d_intermediate.data().get(), d_y.data().get(), s * sizeof(float), cudaMemcpyDeviceToDevice);
        reduce1<<<blocks, threads>>>(d_intermediate.data().get(), d_y.data().get(), s);

        s = (s + threads - 1) / threads;
    }
    cudaEventRecord(stop);

    cudaEventSynchronize(stop);
    milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, start, stop);
    std::cout << 2 * size / 1e9 / (milliseconds / 1000.0) << " GFLOPS, " << milliseconds << " ms" << std::endl;

    res = d_y[0];

    // Verify the result & print first error if any
    correct = (res == size * a);

    std::cout << res << std::endl;
    std::cout << size * a << std::endl;

    if (correct) {
        std::cout << "Reduce test passed!" << std::endl;
    } else {
        std::cout << "Reduce test failed!" << std::endl;
    }

    // cudaDeviceProp prop{};
    // int device;
    // cudaGetDevice(&device);
    // cudaGetDeviceProperties(&prop, device);
    //
    // std::cout << "Device " << prop.name << std::endl;
    // std::cout << prop.major << "." << prop.minor << std::endl;
    // std::cout << prop.maxThreadsPerBlock << std::endl;
    // std::cout << prop.maxGridSize[0] <<std::endl;
    // std::cout << prop.maxBlocksPerMultiProcessor << std::endl;
    // std::cout << prop.maxThreadsPerMultiProcessor << std::endl;


    return 0;
}