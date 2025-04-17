#ifndef KERNELS_CUH
#define KERNELS_CUH

__global__ void axpy(float a, const float *x, const float *y, int size, float *res);

__global__ void reduce0(const float *x, float *y, int N);

__global__ void reduce1(const float *g_idata, float *g_odata, unsigned int n);

__global__ void spMV(const int *ia, const int *ja, const float *a, const float *x, float *y, int size);

__global__ void dot(const float *x, const float *y, float *z, int N);

__global__ void multiply(const float *x, const float *y, float *z, int N);

#endif //KERNELS_CUH
