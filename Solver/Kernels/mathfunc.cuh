#ifndef KERNELS_CUH
#define KERNELS_CUH

__global__ void axpy(double a, double *x, double *y, int size, double *res);

#endif //KERNELS_CUH
