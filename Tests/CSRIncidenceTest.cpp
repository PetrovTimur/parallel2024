#include "Solver/csr.h"
#include "Utilities/input.h"
#include <iostream>
#include <vector>


int main() {
    int Nx = 3;
    int Ny = 2;

    auto stats = input(Nx, Ny, 1, 1); // total_cells, nodes, total_split_cells, total_edges, ...

    // int Ne = Nx * Ny;
    // int Nn = (Nx + 1) * (Ny + 1);
    // int nnz = 4 * Nx * Ny;

    printArray(stats, 5);

    int Ne = stats[0] - stats[2] + 2 * stats[2];
    int Nn = stats[1];
    int nnz = (stats[0] - stats[2]) * 4 + 2 * stats[2] * 3;

    auto csr = makeIncidenceMatrixCSR(Nx, Ny, 1, 1, Ne, nnz);
    auto ia_en = csr.first;
    auto ja_en = csr.second;

    printArray(ia_en, Ne + 1);
    printArray(ja_en, nnz);

    auto csr_transposed = transposeCSR(ia_en, ja_en, Ne, Nn);
    auto ia_ne = csr_transposed.first;
    auto ja_ne = csr_transposed.second;

    auto csr_adjacent = buildAdjacencyMatrixCSR(ia_ne, ja_ne, Ne, Nn);
    auto ia_adj = csr_adjacent.first;
    auto ja_adj = csr_adjacent.second;

    printArray(ia_adj, Ne + 1);
    printArray(ja_adj, ia_adj[Ne]);

    delete[] ia_en;
    delete[] ja_en;
    delete[] ia_ne;
    delete[] ja_ne;
    delete[] ia_adj;
    delete[] ja_adj;

    return 0;
}