#include "mathfunc.h"

double dot(std::vector<double> &x, std::vector<double> &y) {
    double sum = 0;

    #pragma omp parallel for reduction(+:sum)
    for (int i = 0; i < x.size(); i++) {
        sum += x[i] * y[i];
    }

    return sum;
}

void spMV(std::vector<int> &ia, std::vector<int> &ja, std::vector<double> &a, std::vector<double> &b,
    std::vector<double> &res) {
    #pragma omp parallel for
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

void axpy(double a, std::vector<double> x, std::vector<double> y, std::vector<double> &res) {
    #pragma omp parallel for
    for (int i = 0; i < x.size(); i++) {
        res[i] = a * x[i] + y[i];
    }
}
