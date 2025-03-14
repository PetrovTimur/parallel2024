#ifndef INPUT_H
#define INPUT_H
#include <iostream>
#include <tuple>
#include <vector>

#ifdef USE_MPI
std::tuple<int, int, int, int, int, int, int, int, int, int> input(int Nx, int Ny, int Px, int Py, int MyID, std::vector<int> &L2G, std::vector<int> &G2L, std::vector<int> &Part);
#else
int *input(int Nx, int Ny, int K1, int K2);
#endif

std::tuple<int, int, int, int*> readData(const std::string &elementsTxtPath, const std::string &elementsDatPath);

template <typename T>
void printMatrix(std::vector<std::vector<T>> &matrix, std::ostream& os = std::cout) {
    for (int i = 0; i < matrix.size(); i++) {
        for (int j = 0; j < matrix[i].size(); j++) {
            os << matrix[i][j] << " ";
        }
        os << std::endl;
    }
}

template <typename T>
void printVector(std::vector<T> &a, std::ostream& os = std::cout) {
    for (int i = 0; i < a.size(); i++) {
        os << a[i] << " ";
    }
    os << std::endl;
}

template <typename T>
void printArray(T* a, const int size, std::ostream& os = std::cout) {
    for (int i = 0; i < size; i++) {
        os << a[i] << " ";
    }
    os << std::endl;
}

#endif //INPUT_H
