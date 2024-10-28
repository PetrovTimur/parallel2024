#include <iostream>
#include <string>
#include <tuple>
#include <vector>
#include <cmath>
#include <omp.h>
#include "csr.h"
#include "input.h"
#include "solvers.h"

int main(int argc, char** argv) {
    double start = omp_get_wtime();

    omp_set_num_threads(omp_get_max_threads());

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
    std::vector<double> diag(nodes);

    makeCSR(Nx, Ny, K1, K2, ia, ja);
    fillCSR(ia, ja, a, b, diag);

    std::vector<double> res(nodes);

    solve(ia, ja, a, b, diag, res);

    double end = omp_get_wtime();
    printf("Work took %f seconds\n", end - start);

    /*for (int i = 0; i < nodes; i++) {
        std::cout << res[i] << " ";
    }*/

    std::cout << res[253] << std::endl;

    return 0;
}
