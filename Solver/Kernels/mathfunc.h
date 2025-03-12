#ifndef MATHFUNC_H
#define MATHFUNC_H
#include <vector>

void dot(const double *x, const double *y, int size, double &res);

#ifdef USE_MPI
void spMV(std::vector<int> &ia, std::vector<int> &ja, std::vector<double> &a, std::vector<double> &b, std::vector<double> &b_halo, std::vector<double> &res);
#else
void spMV(const int *ia, const int *ja, const double *a, const double *b, int size, double *res);
#endif

void axpy(double a, const double *x, const double *y, int size, double *res);

void arrCopy(double *y, const double *x, int size);

void arrInit(double *x, double a, int size);

#endif //MATHFUNC_H
