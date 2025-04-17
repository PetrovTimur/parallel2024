#ifndef SOLVERS_CUH
#define SOLVERS_CUH

int solve(const int *ia, const int *ja, const float *a, const float *b, const float *diag, int size, float *res, double eps, int maxit);

#endif //SOLVERS_CUH
