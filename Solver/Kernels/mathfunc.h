#ifndef MATHFUNC_H
#define MATHFUNC_H
#include <vector>

void dot(std::vector<double> &x, std::vector<double> &y, double &res);

void spMV(std::vector<int> &ia, std::vector<int> &ja, std::vector<double> &a, std::vector<double> &b, std::vector<double> &res);

void axpy(double a, std::vector<double> x, std::vector<double> y, std::vector<double> &res);

void vecCopy(std::vector<double> &x, std::vector<double> &y);

void vecInit(std::vector<double> &x, double a);

#endif //MATHFUNC_H
