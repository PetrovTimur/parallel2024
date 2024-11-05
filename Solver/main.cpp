#include <iostream>
#include <string>
#include <tuple>
#include <vector>
#include <fstream>
#include <omp.h>
#include "csr.h"
#include "Utilities/input.h"
#include "solvers.h"
#include "Utilities/argparse.h"


int main(int argc, char** argv) {
    struct arguments arguments{};

    /* Default values. */
    arguments.output_file = "-";

    /* Parse our arguments; every option seen by parse_opt will
       be reflected in arguments. */
    argp_parse (&argp, argc, argv, 0, nullptr, &arguments);

    omp_set_num_threads(omp_get_max_threads());

    std::fstream out;
    if (arguments.output_file[0] != '-') {
        out.open(arguments.output_file, std::ios::out);
        std::cout.rdbuf(out.rdbuf()); //save and redirect
    }

    int Nx = std::stoi(arguments.args[0]);
    int Ny = std::stoi(arguments.args[1]);
    int K1 = std::stoi(arguments.args[2]);
    int K2 = std::stoi(arguments.args[3]);

    std::cout << "Nx = " << Nx << std::endl;

    double start = omp_get_wtime();

    std::tuple<int, int> t = input(Nx, Ny, K1, K2);
    int nodes = std::get<0>(t);
    int nonzero_elements =std::get<1>(t);

    std::vector<int> ia(nodes + 1);
    std::vector<int> ja(nonzero_elements);
    std::vector<double> a(nonzero_elements);
    std::vector<double> b(nodes);
    std::vector<double> diag(nodes);

    makeCSR(Nx, Ny, K1, K2, ia, ja);

    #ifdef USE_DEBUG_MODE
    std::vector<std::vector<bool>> matrix(nodes + 1, std::vector<bool>(nodes + 1, false));
    buildAdjacencyMatrix(ia, ja, matrix);
    printMatrix(matrix);
    #endif

    fillCSR(ia, ja, a, b, diag);

    #ifdef USE_DEBUG_MODE
    std::cout << "IA: ";
    printVector(ia);

    std::cout << "JA: ";
    printVector(ja);

    std::cout << "A: ";
    printVector(a);
    #endif
    std::vector<double> res(nodes);

    int iterations = solve(ia, ja, a, b, diag, res);

    double end = omp_get_wtime();
    std::cout << "Work took " << end - start << " seconds\n";
    std::cout << "Convergence required "  << iterations << " iterations\n";

    #ifdef USE_DEBUG_MODE
    std::cout << "res: ";
    printVector(res);
    #endif


    out.close();

    return 0;
}
