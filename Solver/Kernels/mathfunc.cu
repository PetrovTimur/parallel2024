#include "cuda_runtime.h"
#include "thrust/device_vector.h"

__global__ void axpy(double a, double *x, double *y, int size, double *res) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < size) {
        res[i] = a * x[i] + y[i];
    }
}