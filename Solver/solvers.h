#ifndef SOLVERS_H
#define SOLVERS_H

#ifdef USE_MPI
#include <mpi.h>
#endif
#include <vector>

#ifdef USE_MPI
int solve(int MyID, const std::vector<int> &Part, const std::vector<int> &L2G, const std::vector<int> &ia, const std::vector<int> &ja, const std::vector<double> &a, const std::vector<double> &b,
          const std::vector<double> &diag, std::vector<double> &res, double eps, int maxit);
#else
int solve(const int *ia, const int *ja, const double *a, const double *b,
          const double *diag, int size, double *res, double eps, int maxIterations);
#endif


#endif //SOLVERS_H
