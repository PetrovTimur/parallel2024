#include <iostream>
#include <string>
#include <tuple>
#include <vector>
#include <cmath>

std::tuple<int, int> input(int Nx, int Ny, int K1, int K2) {
    std::cout << "Hello, World!" << std::endl;
    std::cout << Nx << " " << Ny << " " << K1 << " " << K2 << std::endl;

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

/*
flat_idx = i * Nx + j;
block_idx = (flat_idx) % (K1 + K2)
if (block_idx < K1)
    return false
else
    return true
*/

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

void fillCSR(std::vector<int> &ia, std::vector<int> &ja, std::vector<double> &a, std::vector<double> &b) {
    for (int i = 0; i < ia.size() - 1; i++) {
        double sum = 0;
        int k_i;
        for (int k = ia[i]; k < ia[i + 1]; k++) {
            int j = ja[k];
            if (i != j) {
                a[k] = std::cos(i * j + i + j);
                sum += a[k];
            } else
                k_i = k;
            std::cout << "el: " << i << ", j: " << j << ", val = " << a[k] << std::endl;
        }

        a[k_i] = 1.234 * sum;
        std::cout << "el: " << i << ", j: " << ja[k_i] << ", val = " << a[k_i] << std::endl;
    }
}

double dot(std::vector<double> &x, std::vector<double> &y) {
    double sum = 0;

    // #pragma omp parallel for reduction(+:sum)
    for (int i = 0; i < x.size(); i++) {
        sum += x[i] * y[i];
    }

    return sum;
}

void solve() {

}

int main(int argc, char** argv) {

    int Nx = std::stoi(argv[1]);
    int Ny = std::stoi(argv[2]);
    int K1 = std::stoi(argv[3]);
    int K2 = std::stoi(argv[4]);

    std::tuple<int, int> t = input(Nx, Ny, K1, K2);
    int nodes = std::get<0>(t);
    int nonzero_elements =std::get<1>(t);

    std::vector<int> ia(nodes + 1);
    std::vector<int> ja(nonzero_elements);
    std::vector<double> a(nonzero_elements);
    std::vector<double> b(nodes);

    makeCSR(Nx, Ny, K1, K2, ia, ja);
    fillCSR(ia, ja, a, b);

    return 0;
}
