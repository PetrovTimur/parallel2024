#include "thrust/device_vector.h"
#include <thrust/host_vector.h>
#include "Solver/Kernels/mathfunc.cuh"
#include "cuda_runtime.h"
#include <iostream>


int main() {
    int blocks = 128;
    int threads = 256;
    int size = 100000000;
    float a = 2.f;
    float b = 3.f;

    thrust::device_vector<float> d_x(size, a);
    thrust::device_vector<float> d_y(size, b);
    thrust::device_vector<float> d_z(blocks);

    float res_host;
    float *res;
    cudaMalloc(&res, sizeof(float));

    dot_gpu(threads, blocks, d_x.data().get(), d_y.data().get(), d_z.data().get(), res, size);

    cudaMemcpy(&res_host, res, sizeof(float), cudaMemcpyDeviceToHost);

    // Verify the result & print first error if any
    bool correct = (res_host == a * b * size);

    std::cout << res_host << std::endl;
    std::cout << size * a * b << std::endl;


    if (correct) {
        std::cout << "Reduce test passed!" << std::endl;
    } else {
        std::cout << "Reduce test failed!" << std::endl;
    }

    cudaFree(res);

    return 0;
}