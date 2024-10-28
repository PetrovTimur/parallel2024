#ifndef INPUT_H
#define INPUT_H
#include <iostream>
#include <tuple>
#include <vector>

std::tuple<int, int> input(int Nx, int Ny, int K1, int K2);

template <typename T>
void printMatrix(std::vector<std::vector<T>> &matrix) {
    for (int i = 0; i < matrix.size(); i++) {
        for (int j = 0; j < matrix[i].size(); j++) {
            std::cout << matrix[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

#endif //INPUT_H
