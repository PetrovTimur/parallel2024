#include "cuda_runtime.h"
#include "thrust/device_vector.h"

__global__ void axpy(const double a, const double *x, const double *y, const int size, double *res) {
    const unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < size) {
        res[i] = a * x[i] + y[i];
    }
}

__global__ void reduce0(const float *x, float *y, const int N) {
    extern __shared__ float tsum[];

    const unsigned int id = threadIdx.x;
    const unsigned int tid = blockIdx.x * blockDim.x + threadIdx.x;
    const unsigned int stride = gridDim.x * blockDim.x;

    tsum[id] = 0.0f;

    for (unsigned int k = tid; k < N; k += stride) {
        tsum[id] += x[k];
    }

    __syncthreads();

    for (unsigned int k = blockDim.x / 2; k > 0; k /= 2) {
        if (id < k) {
            tsum[id] += tsum[id + k];
        }
        __syncthreads();
    }

    if (id == 0) {
        y[blockIdx.x] = tsum[0];
    }
}

__global__ void reduce1(const float *g_idata, float *g_odata, const unsigned int n) {
    extern __shared__ float sdata[];

    // load shared mem
    const unsigned int tid = threadIdx.x;
    const unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

    sdata[tid] = (i < n) ? g_idata[i] : 0;

    __syncthreads();


    // do reduction in shared mem
    for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (tid < s) {
            sdata[tid] += sdata[tid + s];
        }

        __syncthreads();
    }

    // write result for this block to global mem
    if (tid == 0) g_odata[blockIdx.x] = sdata[0];
}
