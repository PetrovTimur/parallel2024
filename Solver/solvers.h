#ifndef SOLVERS_H
#define SOLVERS_H
#ifdef USE_MPI
#include <mpi.h>
#endif
#include <vector>

#ifdef USE_MPI
int solve(int MyID, int Px, int top_halo, int left_halo, int right_halo, int bottom_halo, int i_count, int j_count,
        std::vector<int> &ia, std::vector<int> &ja, std::vector<double> &a, std::vector<double> &b,
        std::vector<double> &diag, std::vector<double> &res);
#else
int solve(const int *ia, const int *ja, const double *a, const double *b,
          const double *diag, int size, double *res, double eps, int maxit);
#endif


#endif //SOLVERS_H
