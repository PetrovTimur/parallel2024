#include <iostream>
#include <string>
#include <tuple>
#include <vector>
#include <cmath>
#include <fstream>
#include <omp.h>
#include "csr.h"
#include "Utilities/input.h"
#include "solvers.h"
#include "Utilities/argparse.h"


int main(int argc, char** argv) {
    struct arguments arguments{};

    /* Default values. */
    arguments.verbose = 0;
    arguments.output_file = "-";

    /* Parse our arguments; every option seen by parse_opt will
       be reflected in arguments. */
    argp_parse (&argp, argc, argv, 0, nullptr, &arguments);

    // printf ("ARG1 = %s\nARG2 = %s\nARG3 = %s\nARG4 = %s\nOUTPUT_FILE = %s\n"
    //         "VERBOSE = %s\n",
    //         arguments.args[0], arguments.args[1],
    //         arguments.args[2], arguments.args[3],
    //         arguments.output_file,
    //         arguments.verbose ? "yes" : "no");

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

    // std::vector<std::vector<bool>> matrix(nodes + 1, std::vector<bool>(nodes + 1, false));

    // buildAM(ia, ja, matrix);

    // printMatrix(matrix);

    fillCSR(ia, ja, a, b, diag);

    std::cout << "IA: ";

    for (auto x: ia) {
        std::cout << x << " ";
    }
    std::cout << std::endl << "JA: ";

    for (auto y : ja) {
        std::cout << y << " ";
    }
    std::cout << std::endl << "A: ";

    for (auto z : a) {
        std::cout << z << " ";
    }
    std::cout << std::endl;

    std::vector<double> res(nodes);

    // Nx = 2, Ny = 2, K1 = 1, K2 = 1
    solve(ia, ja, a, b, diag, res);

    double end = omp_get_wtime();
    std::cout << "Work took " << end - start << " seconds\n";

    for (int i = 0; i < nodes; i++) {
        std::cout << res[i] << " ";
    }

    std::cout << std::endl;

    out.close();

    return 0;
}
