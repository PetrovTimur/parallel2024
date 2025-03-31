#ifndef CSR_H
#define CSR_H
#include <vector>

#ifdef USE_MPI
void makeCSR(int Nx, int Ny, int K1, int K2, int i_start, int i_end, int j_start, int j_end, std::vector<int> &G2L, std::vector<int> &ia, std::vector<int> &ja);
#else
void makeCSR(int Nx, int Ny, int K1, int K2, std::vector<int> &ia, std::vector<int> &ja);
#endif

std::pair<int*, int*> makeIncidenceMatrixCSR(int Nx, int Ny, int K1, int K2, int Ne, int nnz);

#ifdef USE_MPI
void fillCSR(std::vector<int> &ia, std::vector<int> &ja, std::vector<int> &L2G, std::vector<double> &a, std::vector<double> &b,
    std::vector<double> &diag);
#else
void fillCSR(const int *ia, const int *ja, double *a, double *b, double *diag, int size);
#endif

void buildMatrixFromCSR(std::vector<int> &ia, std::vector<int> &ja, std::vector<std::vector<bool>> &matrix);

void buildFilledMatrix(std::vector<int> &ia, std::vector<int> &ja, std::vector<double> &a, std::vector<std::vector<double>> &matrix);

std::pair<int*, int*> transposeCSR(const int* ia, const int* ja, int Ne, int Nn);

std::pair<int*, int*> buildAdjacencyMatrixCSR(const int* ia_ne, const int* ja_ne, int Ne, int Nn);

std::pair<int *, int *> buildAdjacencyMatrixCSRUsingSets(const int *ia_en, const int *ja_en, const int *ia_ne, const int *ja_ne, const int Ne, const int Nn);

std::pair<int*, int*> buildAdjacencyMatrixCSRUsingSort(const int *ia_en, const int *ja_en, const int *ia_ne, const int *ja_ne, const int Ne, const int Nn);


#endif //CSR_H
