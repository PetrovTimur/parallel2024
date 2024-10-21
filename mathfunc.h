#ifndef MATHFUNC_H
#define MATHFUNC_H
#include <vector>

double dot(std::vector<double> &x, std::vector<double> &y);

void spMV(std::vector<int> &ia, std::vector<int> &ja, std::vector<double> &a, std::vector<double> &b, std::vector<double> &res);

#endif //MATHFUNC_H
