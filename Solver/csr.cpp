#include "csr.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <ostream>

#ifdef USE_MPI
void makeCSR(int Nx, int Ny, int K1, int K2, int i_start, int i_end, int j_start, int j_end, std::vector<int> &G2L, std::vector<int> &ia, std::vector<int> &ja) {
    int cycle_length = K1 + K2;
    // std::vector<int> extras(Ny + 1);
    // #pragma omp parallel for proc_bind(spread)
    // for (int i = i_start - (i_start > 0); i <= i_end; i++) {
    //     extras[i] = i * Nx / cycle_length * K2 + std::max(0, i * Nx % cycle_length - K1);
    //     // std::cout << extras[i] << " ";
    // }
    // std::cout << std::endl;

    int k = 0;
    // #pragma omp parallel for proc_bind(spread)
    for (int i = i_start; i <= i_end; i++) {
        // std::cout << i << " " << k << std::endl;
        // int k = 2 * Nx * i + (Nx + 1) * std::max(0, 2 * i - 1) + i * (Nx + 1) + (i > 0 ? extras[i] + extras[i - 1] : 0);
        // std::cout << i << " " << k << std::endl;
        for (int j = j_start; j <= j_end; j++) {
            int flat_idx = i * (Nx + 1) + j;
            // std::cout << i << " " << j << " " << k << " " << flat_idx << " " << G2L[flat_idx] << std::endl;
            // int top_idx = (i - 1) * (Nx + 1) + j;
            // int right_idx = i * (Nx + 1) + j + 1;
            // int bottom_idx = (i + 1) * (Nx + 1) + j;
            // int left_idx = i * (Nx + 1) + j - 1;
            ia[G2L[flat_idx]] = k;

            if (i != 0)
                ja[k++] = G2L[(i - 1) * (Nx + 1) + j];
            if (i != 0 && j != Nx && ((i - 1) * Nx + j) % cycle_length >= K1)
                ja[k++] = G2L[(i - 1) * (Nx + 1) + j + 1];
            if (j != 0)
                ja[k++] = G2L[i * (Nx + 1) + j - 1];

            ja[k++] = G2L[flat_idx];

            if (j != Nx)
                ja[k++] = G2L[i * (Nx + 1) + j + 1];
            if  (j != 0 && i != Ny && (i * Nx + j - 1) % cycle_length >= K1)
                ja[k++] = G2L[(i + 1) * (Nx + 1) + j - 1];
            if (i != Ny)
                ja[k++] = G2L[(i + 1) * (Nx + 1) + j];
        }
    }
    // std::cout << "hiiiii" << std::endl;
    ia[ia.size() - 1] = k;
    // std::cout << "finally" << std::endl;
    ja.resize(k);
}

#else
void makeCSR(int Nx, int Ny, int K1, int K2, std::vector<int> &ia, std::vector<int> &ja) {
    int cycle_length = K1 + K2;
    std::vector<int> extras(Ny + 1);
#pragma omp parallel for proc_bind(spread)
    for (int i = 0; i < Ny + 1; i++) {
        extras[i] = i * Nx / cycle_length * K2 + std::max(0, i * Nx % cycle_length - K1);
        // std::cout << extras[i] << " ";
    }
    // std::cout << std::endl;

    // int k = 0;
#pragma omp parallel for proc_bind(spread)
    for (int i = 0; i < Ny + 1; i++) {
        // std::cout << i << " " << k << std::endl;
        int k = 2 * Nx * i + (Nx + 1) * std::max(0, 2 * i - 1) + i * (Nx + 1) + (i > 0 ? extras[i] + extras[i - 1] : 0);
        // std::cout << i << " " << k << std::endl;
        for (int j = 0; j < Nx + 1; j++) {
            int flat_idx = i * (Nx + 1) + j;
            // int top_idx = (i - 1) * (Nx + 1) + j;
            // int right_idx = i * (Nx + 1) + j + 1;
            // int bottom_idx = (i + 1) * (Nx + 1) + j;
            // int left_idx = i * (Nx + 1) + j - 1;
            ia[flat_idx] = k;

            if (i != 0)
                ja[k++] = (i - 1) * (Nx + 1) + j;
            if (i != 0 && j != Nx && ((i - 1) * Nx + j) % cycle_length >= K1)
                ja[k++] = (i - 1) * (Nx + 1) + j + 1;
            if (j != 0)
                ja[k++] = i * (Nx + 1) + j - 1;

            ja[k++] = flat_idx;

            if (j != Nx)
                ja[k++] = i * (Nx + 1) + j + 1;
            if  (j != 0 && i != Ny && (i * Nx + j - 1) % cycle_length >= K1)
                ja[k++] = (i + 1) * (Nx + 1) + j - 1;
            if (i != Ny)
                ja[k++] = (i + 1) * (Nx + 1) + j;
        }
    }
    ia[(Nx + 1) * (Ny + 1)] = ja.size();
}

#endif


#ifdef USE_MPI
void fillCSR(std::vector<int> &ia, std::vector<int> &ja, std::vector<int> &L2G, std::vector<double> &a, std::vector<double> &b,
    std::vector<double> &diag) {
    // #pragma omp parallel for proc_bind(spread)
    for (int i = 0; i < ia.size() - 1; i++) {
        double sum = 0;
        int k_i = 0;
        for (int k = ia[i]; k < ia[i + 1]; k++) {
            int j = ja[k];
            if (i != j) {
                a[k] = std::cos(L2G[i] * L2G[j] + L2G[i] + L2G[j]);
                sum += std::abs(a[k]);
            } else
                k_i = k;
            // std::cout << "el: " << i << ", j: " << j << ", val = " << a[k] << std::endl;
        }

        a[k_i] = 1.234 * sum;
        diag[i] = a[k_i];
        // std::cout << "el: " << i << ", j: " << ja[k_i] << ", val = " << a[k_i] << std::endl;
    }

    // #pragma omp parallel for proc_bind(spread)
    for (int i = 0; i < b.size(); i++)
        b[i] = std::sin(L2G[i]);
}
#else
void fillCSR(std::vector<int> &ia, std::vector<int> &ja, std::vector<double> &a, std::vector<double> &b,
    std::vector<double> &diag) {
#pragma omp parallel for proc_bind(spread)
    for (int i = 0; i < ia.size() - 1; i++) {
        double sum = 0;
        int k_i = 0;
        for (int k = ia[i]; k < ia[i + 1]; k++) {
            int j = ja[k];
            if (i != j) {
                a[k] = std::cos(i * j + i + j);
                sum += std::abs(a[k]);
            } else
                k_i = k;
            // std::cout << "el: " << i << ", j: " << j << ", val = " << a[k] << std::endl;
        }

        a[k_i] = 1.234 * sum;
        diag[i] = a[k_i];
        // std::cout << "el: " << i << ", j: " << ja[k_i] << ", val = " << a[k_i] << std::endl;
    }

#pragma omp parallel for proc_bind(spread)
    for (int i = 0; i < b.size(); i++)
        b[i] = std::sin(i);
}
#endif

void buildAdjacencyMatrix(std::vector<int> &ia, std::vector<int> &ja, std::vector<std::vector<bool>> &matrix) {
    for (int i = 0; i < ia.size() - 1; i++) {
        for (int k = ia[i]; k < ia[i + 1]; k++) {
            matrix[i][ja[k]] = true;
        }
    }
}

void buildMatrix(std::vector<int> &ia, std::vector<int> &ja, std::vector<double> &a,
    std::vector<std::vector<double>> &matrix) {
    for (int i = 0; i < ia.size() - 1; i++) {
        for (int k = ia[i]; k < ia[i + 1]; k++) {
            matrix[i][ja[k]] = a[k];
        }
    }
}
