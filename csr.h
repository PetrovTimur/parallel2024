#ifndef CSR_H
#define CSR_H
#include <vector>

void makeCSR(int Nx, int Ny, int K1, int K2, std::vector<int> &ia, std::vector<int> &ja);

void fillCSR(std::vector<int> &ia, std::vector<int> &ja, std::vector<double> &a, std::vector<double> &b, std::vector<double> &diag);

void buildAdjacencyMatrix(std::vector<int> &ia, std::vector<int> &ja, std::vector<std::vector<bool>> &matrix);

void buildMatrix(std::vector<int> &ia, std::vector<int> &ja, std::vector<double> &a, std::vector<std::vector<double>> &matrix);

#endif //CSR_H
