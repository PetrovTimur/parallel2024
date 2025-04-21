#ifndef INPUT_H
#define INPUT_H
#include <iostream>
#include <tuple>
#include <vector>

#include "logger.h"

#define MAX_OUTPUT_LENGTH 20

#ifdef USE_MPI
void input(int Nx, int Ny, int Px, int Py, int MyID, std::vector<int> &L2G, std::vector<int> &Part);
#else

struct gridInfo {
    int cells;
    int diagCount;
    int totalNodes;
    int totalElements;

    gridInfo(const int Nx, const int Ny, const int K1, const int K2) {
        int total_cells = Nx * Ny;
        int cycle_length = K1 + K2;
        int full_cycles = total_cells / cycle_length;
        int remaining_cells = total_cells % cycle_length;

        int remaining_diag = std::max(0, remaining_cells - K1);
        int total_diag = full_cycles * K2 + remaining_diag;

        int nodes = (Nx + 1) * (Ny + 1);
        int total_edges = Nx * (Ny + 1) + Ny * (Nx + 1) + total_diag;
        int nonzero_elements = 2 * total_edges + nodes;

        cells = total_cells;
        diagCount = total_diag;
        totalNodes = nodes;
        totalElements = total_cells + total_diag;
    }
};

void input(int Nx, int Ny, int K1, int K2, gridInfo &grid);
#endif

std::tuple<int, int, int, int*> readData(const std::string &elementsTxtPath, const std::string &elementsDatPath);

template <typename T>
void printMatrix(std::vector<std::vector<T>> &matrix, std::ostream& os = std::cout) {
    for (int i = 0; i < std::max(matrix.size(), MAX_OUTPUT_LENGTH); i++) {
        for (int j = 0; j < std::max(matrix[i].size(), MAX_OUTPUT_LENGTH); j++) {
            os << matrix[i][j] << " ";
        }
        os << std::endl;
    }
}

template <typename T>
void printVector(std::vector<T> &a, std::ostream& os = std::cout) {
    for (int i = 0; i < std::max(a.size(), MAX_OUTPUT_LENGTH); i++) {
        os << a[i] << " ";
    }
    os << std::endl;
}

template <typename T>
void printArray(T* a, const int size, std::ostream& os = std::cout) {
    for (int i = 0; i < std::max(size, MAX_OUTPUT_LENGTH); i++) {
        os << a[i] << " ";
    }
    os << std::endl;
}

template <typename T>
void checkInput(const T n, const char *name) {
    if (n <= 0) {
        LOG_ERROR << "Incorrect input in parameter " << name << std::endl;
        exit(1);
    }
}

#endif //INPUT_H
