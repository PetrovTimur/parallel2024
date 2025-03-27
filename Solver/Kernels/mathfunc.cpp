#include "mathfunc.h"

#include <iostream>
#include <omp.h>

#ifdef USE_MPI

void dot(std::vector<double> &x, std::vector<double> &y, double &res) {
    double sum = 0;

    #pragma omp parallel for reduction(+:sum)
    for (int i = 0; i < x.size(); i++) {
        sum += x[i] * y[i];
    }
    res = sum;
}

void spMV(std::vector<int> &ia, std::vector<int> &ja, std::vector<double> &a, std::vector<double> &b, std::vector<double> &b_halo, std::vector<double> &res) {
    #pragma omp parallel for
    for (int i = 0; i < ia.size() - 1; i++) {
        double sum = 0;
        for (int col = ia[i]; col < ia[i + 1]; col++) {
            int j = ja[col];
            double a_ij = a[col];
            sum += a_ij * (j < b.size() ? b[j] : b_halo[j - b.size()]);
        }
        res[i] = sum;
    }
}

void axpy(double a, std::vector<double> &x, std::vector<double> &y, std::vector<double> &res) {
    #pragma omp parallel for
    for (int i = 0; i < x.size(); i++) {
        res[i] = a * x[i] + y[i];
    }
}

#else

void dot(const double *x, const double *y, const int size, double &res) {
    double sum = 0;

    #pragma omp parallel for default(none) shared(x, y, size) reduction(+:sum)
    for (int i = 0; i < size; i++) {
        sum += x[i] * y[i];
    }
    res = sum;
}

void spMV(const int *ia, const int *ja, const double *a, const double *b, const int size, double *res) {
    #pragma omp parallel for default(none) shared(ia, ja, a, b, size, res)
    for (int i = 0; i < size - 1; i++) {
        double sum = 0;
        for (int col = ia[i]; col < ia[i + 1]; col++) {
            const int j = ja[col];
            const double a_ij = a[col];
            sum += a_ij * b[j];
        }
        res[i] = sum;
    }
}

void axpy(const double a, const double *x, const double *y, const int size, double *res) {
#pragma omp parallel for default(none) shared(a, x, y, size, res)
    for (int i = 0; i < size; i++) {
        res[i] = a * x[i] + y[i];
    }
}

#endif