#include "input.h"
#include <fstream>

#ifndef USE_MPI
std::tuple<int, int, int, int *> readData(const std::string &elementsTxtPath, const std::string &elementsDatPath) {
    // Read the first file
    std::ifstream elementsTxt(elementsTxtPath);
    if (!elementsTxt.is_open()) {
        throw std::runtime_error("Failed to open elements.txt");
    }

    int Nn, Ne, M;
    std::string dummy;
    elementsTxt >> dummy >> Nn;
    elementsTxt >> dummy >> Ne;
    elementsTxt >> dummy >> M;
    elementsTxt.close();

    // Read the second file
    std::ifstream elementsDat(elementsDatPath, std::ios::binary);
    if (!elementsDat.is_open()) {
        throw std::runtime_error("Failed to open elements.dat");
    }

    auto data = new int[Ne * M];
    elementsDat.read(reinterpret_cast<char*>(data), Ne * M * sizeof(int));
    elementsDat.close();

    return std::make_tuple(Nn, Ne, M, data);
}
#else

void input(const int Nx, const int Ny, const int Px, const int Py, const int MyID, std::vector<int> &L2G, std::vector<int> &Part) {
    const int MyID_j = MyID % Px;
    const int MyID_i = MyID / Px;

    int j_start = MyID_j * (Nx / Px) + std::min(Nx % Px, MyID_j);
    int j_count = Nx / Px + (Nx % Px > MyID_j);
    int j_end = j_start + j_count - 1;
    int left_halo = MyID_j > 0; // Halo
    j_start -= left_halo;
    int right_halo = MyID_j < Px - 1; // Halo
    j_end += right_halo;

    int i_start = MyID_i * (Ny / Py) + std::min(Ny % Py, MyID_i);
    int i_count = Ny / Py + (Ny % Py > MyID_i);
    int i_end = i_start + i_count - 1;
    int top_halo = MyID_i > 0; // Halo
    i_start -= top_halo;
    int bottom_halo = MyID_i < Py - 1; // Halo
    i_end += bottom_halo;

    // std::cout << "MyID: " << MyID << std::endl;
    // std::cout << "i_start: " << i_start << ", i_end: " << i_end << std::endl;
    // std::cout << "j_start: " << j_start << ", j_end: " << j_end << std::endl;

    // start/end are for local elements

    // int k = 0;
    // L2G.resize((i_end - i_start + 1) * (j_end - j_start + 1));
    for (int i = i_start + top_halo; i <= i_end - bottom_halo; i++) {
        for (int j = j_start + left_halo; j <= j_end - right_halo; j++) {
            L2G.push_back(i * Nx + j);
            Part.push_back(MyID);
            // L2G[k++] = i * (Nx + 1) + j;
        }
    }

    // int N_local = L2G.size();
    if (top_halo) {
        if (left_halo) {
            L2G.push_back(i_start * (Nx) + j_start);
            // L2G[k++] = i_start * (Nx + 1) + j_start;
            Part.push_back(MyID - Px - 1);
        }

        for (int j = j_start + left_halo; j <= j_end - right_halo; j++) {
            L2G.push_back(i_start * (Nx) + j);
            Part.push_back(MyID - Px);
            // L2G[k++] = i_start * (Nx + 1) + j;
        }
        if (right_halo) {
            L2G.push_back(i_start * (Nx) + j_end);
            // L2G[k++] = i_start * (Nx + 1) + j_end;
            Part.push_back(MyID - Px + 1);
        }
    }
    if (left_halo) {
        for (int i = i_start + top_halo; i <= i_end - bottom_halo; i++) {
            // L2G[k++] = i * (Nx + 1) + j_start;
            L2G.push_back(i * (Nx) + j_start);
            Part.push_back(MyID - 1);
        }
    }
    if (right_halo) {
        for (int i = i_start + top_halo; i <= i_end - bottom_halo; i++) {
            // L2G[k++] = i * (Nx + 1) + j_end;
            L2G.push_back(i * (Nx) + j_end);
            Part.push_back(MyID + 1);
        }
    }
    if (bottom_halo) {
        if (left_halo) {
            // L2G[k++] = i_end * (Nx + 1) + j_start;
            L2G.push_back(i_end * (Nx) + j_start);
            Part.push_back(MyID + Px - 1);
        }

        for (int j = j_start + left_halo; j <= j_end - right_halo; j++) {
            // L2G[k++] = i_end * (Nx + 1) + j;
            L2G.push_back(i_end * (Nx) + j);
            Part.push_back(MyID + Px);
        }

        if (right_halo) {
            // L2G[k++] = i_end * (Nx + 1) + j_end;
            L2G.push_back(i_end * (Nx) + j_end);
            Part.push_back(MyID + Px + 1);
        }
    }

    // L2G.resize(k);
}
#endif
