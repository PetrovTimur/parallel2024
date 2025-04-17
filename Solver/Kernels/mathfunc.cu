#include "cuda_runtime.h"
#include "thrust/device_vector.h"

__global__ void axpy(const float a, const float *x, const float *y, const int size, float *res) {
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

__global__ void spMV(const int *ia, const int *ja, const float *a, const float *x, float *y, const int size) {
    const unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < size) {
        float sum = 0.0;
        for (unsigned int k = ia[i]; k < ia[i + 1]; k++) {
            const unsigned int j = ja[k];
            const float a_ij = a[k];
            sum += x[j] * a_ij;
        }

        y[i] = sum;
    }
}

__global__ void dot(const float *x, const float *y, float *z, const int N) {
    extern __shared__ float tsum[256];

    const unsigned int id = threadIdx.x;
    const unsigned int tid = blockIdx.x * blockDim.x + threadIdx.x;
    const unsigned int stride = gridDim.x * blockDim.x;

    tsum[id] = 0.0f;

    for (unsigned int k = tid; k < N; k += stride) {
        tsum[id] += x[k] * y[k];
    }

    __syncthreads();

    for (unsigned int k = blockDim.x / 2; k > 0; k /= 2) {
        if (id < k) {
            tsum[id] += tsum[id + k];
        }
        __syncthreads();
    }

    if (id == 0) {
        z[blockIdx.x] = tsum[0];
    }
}


__global__ void multiply(const float *x, const float *y, float *z, const int N) {
    const unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    const unsigned int stride = gridDim.x * blockDim.x;

    for (unsigned int k = i; k < N; k += stride) {
        z[k] = x[k] * y[k];
    }
}

