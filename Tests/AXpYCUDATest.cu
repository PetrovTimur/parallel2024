#include "cuda_runtime.h"
#include "thrust/device_vector.h"
#include "Solver/Kernels/mathfunc.cuh"
#include <iostream>
#include <vector>


void axpy_gpu(float a, float *x, float *y, int size, float *res) {
    int threadsPerBlock = 256;
    int blocksPerGrid = (size + threadsPerBlock - 1) / threadsPerBlock;

    float *d_x, *d_y, *d_res;
    cudaMalloc((void**)&d_x, size * sizeof(float));
    cudaMalloc((void**)&d_y, size * sizeof(float));
    cudaMalloc((void**)&d_res, size * sizeof(float));

    cudaMemcpy(d_x, x, size * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_y, y, size * sizeof(float), cudaMemcpyHostToDevice);

    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    cudaEventRecord(start);
    axpy<<<blocksPerGrid, threadsPerBlock>>>(a, d_x, d_y, size, d_res);
    cudaEventRecord(stop);

    cudaEventSynchronize(stop);
    float milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, start, stop);
    std::cout << 2 * size / 1e9 / (milliseconds / 1000.0) << ", " << milliseconds / 1000.0 << std::endl;

    cudaEventDestroy(start);
    cudaEventDestroy(stop);

    cudaMemcpy(res, d_res, size * sizeof(float), cudaMemcpyDeviceToHost);

    cudaFree(d_x);
    cudaFree(d_y);
    cudaFree(d_res);
}

int main() {
    int size = 100000000;
    float a = 2.0;
    std::vector<float> x(size, 1.0);
    std::vector<float> y(size, 2.0);
    std::vector<float> res(size);

    axpy_gpu(a, x.data(), y.data(), size, res.data());

    // Verify the result & print first error if any
    bool correct = true;
    for (int i = 0; i < size; ++i) {
        float expected = a * x[i] + y[i];
        if (res[i] != expected) {
            std::cerr << "Mismatch at index " << i << ": " << res[i] << " != " << expected << std::endl;
            correct = false;
            break;
        }
    }

    if (correct) {
        std::cout << "AXPY test passed!" << std::endl;
    } else {
        std::cout << "AXPY test failed!" << std::endl;
    }

    return 0;
}