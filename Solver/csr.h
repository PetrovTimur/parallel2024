#ifndef CSR_H
#define CSR_H
#include <unordered_map>
#include <vector>

#ifdef USE_MPI
void makeCSR(int Nx, int Ny, int K1, int K2, int i_start, int i_end, int j_start, int j_end, std::vector<int> &G2L, std::vector<int> &ia, std::vector<int> &ja);

void localizeCSR(const int *ia, const int size, int *ja, std::unordered_map<int, int> &G2L);

void constructG2L(std::vector<int> &ia_en, std::vector<int> &ja_en, std::unordered_map<int, int> &G2L);

void fillG2L(std::vector<int> &L2G, std::unordered_map<int, int> &G2L);

void fillCSR(std::vector<int> &ia, std::vector<int> &ja, std::vector<int> &L2G, std::vector<double> &a, std::vector<double> &b,
    std::vector<double> &diag);

void makeIncidenceMatrixCSR(const int Nx, int Ny, const int K1, const int K2, std::vector<int> &L2G, std::vector<int> &ia, std::vector<int> &ja, std
                            ::vector<int> &Part);
#else
void makeCSR(int Nx, int Ny, int K1, int K2, std::vector<int> &ia, std::vector<int> &ja);

template <typename T>
void fillCSR(const int *ia, const int *ja, T *a, T *b, T *diag, int size);

void makeIncidenceMatrixCSR(int Nx, int Ny, int K1, int K2, int i_start, int i_end, int j_start, int j_end,
                            std::vector<int> &ia, std::vector<int>
                            &ja);
#endif

void buildMatrixFromCSR(std::vector<int> &ia, std::vector<int> &ja, std::vector<std::vector<bool>> &matrix);

void buildFilledMatrix(std::vector<int> &ia, std::vector<int> &ja, std::vector<double> &a, std::vector<std::vector<double>> &matrix);

void transposeCSR(std::vector<int> &ia, std::vector<int> &ja, int Nn, int *&ia_new,
                  int *&ja_new);

std::pair<int*, int*> buildAdjacencyMatrixCSR(const int* ia_ne, const int* ja_ne, int Ne, int Nn);

std::pair<int *, int *> buildAdjacencyMatrixCSRUsingSets(const int *ia_en, const int *ja_en, const int *ia_ne, const int *ja_ne, int Ne, int Nn);

#ifdef USE_MPI
void buildAdjacencyMatrixCSRUsingSort(const int *ia_en, const int *ja_en, const int *ia_ne, const int *ja_ne, std::vector<int> &ia_adj, std::vector<int> &
                                      ja_adj, std::vector<int> &Part, int MyID);
#else
std::pair<int*, int*> buildAdjacencyMatrixCSRUsingSort(const int *ia_en, const int *ja_en, const int *ia_ne, const int *ja_ne, int Ne, int Nn);
#endif


#endif //CSR_H
