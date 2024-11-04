#include <omp.h>
#include <iostream>
#include <ostream>
#include <tuple>

#include "Solver/csr.h"
#include "Solver/solvers.h"
#include "Utilities/input.h"

int main() {
    omp_set_num_threads(omp_get_max_threads());
    std::vector<int> ia;
    std::vector<int> ja;
    std::vector<double> a;
    std::vector<double> b;
    std::vector<double> diag;
    std::vector<double> res;

    std::cout << "T = " << omp_get_max_threads() << std::endl;

    for (int k = 1e5; k <= 1e8; k *= 10) {
        int K1 = 30;
        int K2 = 23;
        int Nx = 2000;
        int Ny = k / Nx / 5;
        std::tuple<int, int> t = input(Nx, Ny, K1, K2);
        int nodes = std::get<0>(t);
        int nonzero_elements =std::get<1>(t);

        ia.resize(nodes + 1);
        ja.resize(nonzero_elements);
        a.resize(nonzero_elements);
        b.resize(nodes);
        diag.resize(nodes);
        res.resize(nodes);

        makeCSR(Nx, Ny, K1, K2, ia, ja);
        fillCSR(ia, ja, a, b, diag);

        int runs = 1;
        int iterations = 0;
        double aggregate_time = 0;
        for (int p = 0; p < runs; ++p) {
            // Calculate
            double start = omp_get_wtime();
            iterations = solve(ia, ja, a, b, diag, res);
            double end = omp_get_wtime();

            aggregate_time += end - start;
        }

        double average_time = aggregate_time / runs;

        std::cout << iterations * (3 * nodes + 3 * nodes + a.size()) / (average_time * 1e9) << ", ";
    }
    std::cout << std::endl;
}
