#ifndef SOLVERS_CUH
#define SOLVERS_CUH

void dot_gpu(int threads, int blocks, const float *x, const float *y, float *vec_buf, float *z, int N);
int solve(const int *ia, const int *ja, const float *a, const float *b, const float *diag, int size, float *res, double eps, int maxit);

#endif //SOLVERS_CUH
