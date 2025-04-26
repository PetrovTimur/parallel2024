#include "cuda_runtime.h"
#include "thrust/device_vector.h"
#include "Solver/Kernels/mathfunc.cuh"
#include <iostream>
#include <vector>


void axpy_gpu(float a, float *x, float *y, int size, float *res) {
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    cudaEventRecord(start);
    axpy<<<204, 256>>>(a, x, y, size, res);
    cudaEventRecord(stop);

    cudaEventSynchronize(stop);
    float milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, start, stop);
    std::cout << 2 * size / 1e9 / (milliseconds / 1000.0) << ", " << milliseconds / 1000.0 << std::endl;

    cudaEventDestroy(start);
    cudaEventDestroy(stop);
}

int main(int argc, char *argv[]) {
    int size = std::atoi(argv[1]);
    float a = 2.0;
    std::vector<float> x(size, 1.0);
    std::vector<float> y(size, 2.0);
    std::vector<float> res(size);

    thrust::device_vector<float> d_x = x;
    thrust::device_vector<float> d_y = y;
    thrust::device_vector<float> d_res(size);

    axpy_gpu(a, d_x.data().get(), d_y.data().get(), size, d_res.data().get());

    cudaMemcpy(res.data(), d_res.data().get(), size*sizeof(float), cudaMemcpyDeviceToHost);


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