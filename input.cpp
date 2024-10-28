#include "input.h"
#include <iostream>

std::tuple<int, int> input(int Nx, int Ny, int K1, int K2) {
    // std::cout << "Hello, World!" << std::endl;
    // std::cout << Nx << " " << Ny << " " << K1 << " " << K2 << std::endl;

    int total_cells = Nx * Ny;
    int cycle_length = K1 + K2;
    int full_cycles = total_cells / cycle_length;
    int remaining_cells = total_cells % cycle_length;

    int remaining_ones = std::max(0, remaining_cells - K1);
    int total_ones = full_cycles * K2 + remaining_ones;

    int nodes = (Nx + 1) * (Ny + 1);
    int total_edges = Nx * (Ny + 1) + Ny * (Nx + 1) + total_ones;
    int nonzero_elements = 2 * total_edges + nodes;

    return std::make_tuple(nodes, nonzero_elements);
}

void printMatrix(std::vector<std::vector<bool>> &matrix) {
    for (int i = 0; i < matrix.size(); i++) {
        for (int j = 0; j < matrix[i].size(); j++) {
            std::cout << matrix[i][j] << " ";
        }
        std::cout << std::endl;
    }
}
