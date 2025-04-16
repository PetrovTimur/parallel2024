#ifndef KERNELS_CUH
#define KERNELS_CUH

__global__ void axpy(double a, const double *x, const double *y, int size, double *res);

__global__ void reduce0(const float *x, float *y, int N);

__global__ void reduce1(const float *g_idata, float *g_odata, unsigned int n);

__global__ void spMV(const int *ia, const int *ja, const double *a, const double *x, double *y, int size);

#endif //KERNELS_CUH
