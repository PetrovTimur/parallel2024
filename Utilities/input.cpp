#include "input.h"

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


std::tuple<int, int, int, int, int, int, int, int, int, int> input(int Nx, int Ny, int Px, int Py, int MyID, std::vector<int> &L2G, std::vector<int> &G2L, std::vector<int> &Part) {
    int MyID_j = MyID % Px;
    int MyID_i = MyID / Px;

    int i_start, i_end, i_count, j_start, j_end, j_count;
    j_start = MyID_j * ((Nx + 1) / Px) + std::min((Nx + 1) % Px, MyID_j);
    j_count = (Nx + 1) / Px + ((Nx + 1) % Px > MyID_j);
    j_end = j_start + j_count - 1;
    int left_halo = MyID_j > 0; // Halo
    j_start -= left_halo;
    int right_halo = MyID_j < Px - 1; // Halo
    j_end += right_halo;

    i_start = MyID_i * ((Ny + 1) / Py) + std::min((Ny + 1) % Py, MyID_i);
    i_count = (Ny + 1) / Py + ((Ny + 1) % Py > MyID_i);
    i_end = i_start + i_count - 1;
    int top_halo = MyID_i > 0; // Halo
    i_start -= top_halo;
    int bottom_halo = MyID_i < Py - 1; // Halo
    i_end += bottom_halo;

    // std::cout << "MyID: " << MyID << std::endl;
    // std::cout << "i_start: " << i_start << ", i_end: " << i_end << std::endl;
    // std::cout << "j_start: " << j_start << ", j_end: " << j_end << std::endl;

    int k = 0, N0;
    // std::vector<int> L2G((i_end - i_start + 1) * (j_end - j_start + 1));
    L2G.resize((i_end - i_start + 1) * (j_end - j_start + 1));
    for (int i = i_start + top_halo; i <= i_end - bottom_halo; i++) {
        for (int j = j_start + left_halo; j <= j_end - right_halo; j++) {
            L2G[k++] = i * (Nx + 1) + j;
        }
    }
    N0 = k;
    if (top_halo) {
        for (int j = j_start + left_halo; j <= j_end - right_halo; j++) {
            L2G[k++] = i_start * (Nx + 1) + j;
        }
        if (right_halo)
            L2G[k++] = i_start * (Nx + 1) + j_end;
    }
    if (left_halo) {
        for (int i = i_start + top_halo; i <= i_end - bottom_halo; i++) {
            L2G[k++] = i * (Nx + 1) + j_start;
        }
    }
    if (right_halo) {
        for (int i = i_start + top_halo; i <= i_end - bottom_halo; i++) {
            L2G[k++] = i * (Nx + 1) + j_end;
        }
    }
    if (bottom_halo) {
        if (left_halo)
            L2G[k++] = i_end * (Nx + 1) + j_start;
        for (int j = j_start + left_halo; j <= j_end - right_halo; j++) {
            L2G[k++] = i_end * (Nx + 1) + j;
        }
    }

    L2G.resize(k);

    // if (MyID == 7) {
    //     std::cout << i_start << i_end << j_start << j_end << std::endl;
    //     // std::cout << top_halo << bottom_halo << left_halo << right_halo << std::endl;
    //     std::cout << "L2G size: " << L2G.size() << std::endl;
    //     for (int i = 0; i < L2G.size(); i++) {
    //         std::cout << L2G[i] << " ";
    //     }
    //     std::cout << std::endl;
    // }

    for (int iL = 0; iL < L2G.size(); iL++) {
        G2L[L2G[iL]] = iL;
    }


    for (int i = 0; i < Py; i++) {
        int is = i * ((Ny + 1) / Py) + std::min((Ny + 1) % Py, i);
        int ic = (Ny + 1) / Py + ((Ny + 1) % Py > i) - 1;
        int ie = is + ic;
        for (int j = 0; j < Px; j++) {
            int js = j * ((Nx + 1) / Px) + std::min((Nx + 1) % Px, j);
            int jc = (Nx + 1) / Px + ((Nx + 1) % Px > j) - 1;
            int je = js + jc;

            #ifdef  USE_DEBUG_MODE
            if (MyID == 0) {
                std::cout << "id: " << i * Px + j << std::endl;
                std::cout << is << " " << ie << " " << js << " " << jc << std::endl;
            }
            #endif

            for (int ii = is; ii <= ie; ii++) {
                for (int jj = js; jj <= je; jj++) {
                    Part[ii * (Nx + 1) + jj] = i * Px + j;
                }
            }
        }
    }

    #ifdef USE_DEBUG_MODE
    if (MyID == 0) {
        std::cout << "Part size: " << Part.size() << std::endl;
        for (int i : Part) {
            std::cout << i << " ";
        }
        std::cout << std::endl;
    }
    #endif

    return std::make_tuple(i_start,  i_end, i_count,j_start, j_end, j_count, top_halo, right_halo, bottom_halo, left_halo);
}
