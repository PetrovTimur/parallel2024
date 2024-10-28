#ifndef INPUT_H
#define INPUT_H
#include <tuple>
#include <vector>

std::tuple<int, int> input(int Nx, int Ny, int K1, int K2);

void printMatrix(std::vector<std::vector<bool>> &matrix);

#endif //INPUT_H
