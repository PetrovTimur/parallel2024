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
void dot(std::vector<double> &x, std::vector<double> &y, double &res) {
    double sum = 0;

    #pragma omp parallel for reduction(+:sum) proc_bind(spread)
    for (int i = 0; i < x.size(); i++) {
        sum += x[i] * y[i];
    }
    res = sum;
}

void spMV(std::vector<int> &ia, std::vector<int> &ja, std::vector<double> &a, std::vector<double> &b, std::vector<double> &res) {
    #pragma omp parallel for proc_bind(spread)
    for (int i = 0; i < ia.size() - 1; i++) {
        double sum = 0;
        for (int col = ia[i]; col < ia[i + 1]; col++) {
            int j = ja[col];
            double a_ij = a[col];
            sum += a_ij * b[j];
        }
        res[i] = sum;
    }
}

void axpy(double a, std::vector<double> &x, std::vector<double> &y, std::vector<double> &res) {
    #pragma omp parallel for proc_bind(spread)
    for (int i = 0; i < x.size(); i++) {
        res[i] = a * x[i] + y[i];
    }
}

#endif

void vecCopy(std::vector<double> &x, std::vector<double> &y) {
    #pragma omp parallel for proc_bind(spread)
    for (int i = 0; i < x.size(); i++) {
        y[i] = x[i];
    }
}

void vecInit(std::vector<double> &x, const double a) {
    #pragma omp parallel for proc_bind(spread)
    for (int i = 0; i < x.size(); i++) {
        x[i] = a;
    }
}
