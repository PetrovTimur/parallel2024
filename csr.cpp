#include "csr.h"
#include <cmath>

void makeCSR(int Nx, int Ny, int K1, int K2, std::vector<int> &ia, std::vector<int> &ja) {
    int cycle_length = K1 + K2;
    int k = 0;
    for (int i = 0; i < Ny + 1; i++) {
        for (int j = 0; j < Nx + 1; j++) {
            int flat_idx = i * (Nx + 1) + j;
            // add top, left, bot, right ...
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
    ia[(Nx + 1) * (Ny + 1)] = k;
}

void fillCSR(std::vector<int> &ia, std::vector<int> &ja, std::vector<double> &a, std::vector<double> &b,
    std::vector<double> &diag) {
#pragma omp parallel for
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

    for (int i = 0; i < b.size(); i++)
        b[i] = std::sin(i);
}
