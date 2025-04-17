#include "thrust/device_vector.h"
#include <thrust/host_vector.h>
#include "Solver/Kernels/mathfunc.cuh"
#include "cuda_runtime.h"
#include <iostream>


int main() {
    int blocks = 128;
    int threads = 256;
    int size = 100000000;
    float a = 2.0f;

    thrust::device_vector<float> d_x(size, a);
    thrust::device_vector<float> d_y(size, a);
    thrust::device_vector<float> d_z(blocks, a);

    float res_host;
    float *res;
    cudaMalloc(&res, sizeof(float));

    // dot<<<blocks, threads, threads*sizeof(float)>>>(d_x.data().get(), d_y.data().get(), d_z.data().get(), size);
    // reduce0<<<1, blocks, blocks*sizeof(float)>>>(d_z.data().get(), res, blocks);

    dot_gpu(threads, blocks, d_x.data().get(), d_y.data().get(), d_z.data().get(), res, size);

    cudaMemcpy(&res_host, res, sizeof(float), cudaMemcpyDeviceToHost);

    // Verify the result & print first error if any
    bool correct = (res_host == a * size * a);

    std::cout << res_host << std::endl;
    std::cout << 2 * size * a << std::endl;


    if (correct) {
        std::cout << "Reduce test passed!" << std::endl;
    } else {
        std::cout << "Reduce test failed!" << std::endl;
    }

    cudaFree(res);

    return 0;
}