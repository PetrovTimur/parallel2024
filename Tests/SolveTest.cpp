#include <iostream>
#include <ostream>

#include "solvers.h"
#include "Ops/mathfunc.h"
#include "Utilities/input.h"
// #include "csr.h"


int main() {
    std::vector<int> ia = {0, 5, 10, 15, 20, 25};
    std::vector<int> ja = {
        0, 1, 2, 3, 4,
        0, 1, 2, 3, 4,
        0, 1, 2, 3, 4,
        0, 1, 2, 3, 4,
        0, 1, 2, 3, 4};
    std::vector<double> a = {
        11, 1, 2, 3, 4,
        1, 12, 1, 2, 3,
        2, 1, 13, 1, 2,
        3, 2, 1, 14, 1,
        4, 3, 2, 1, 15};

    std::vector<double> b = {51, 51, 57, 71, 95};
    std::vector<double> diag = {11, 12, 13, 14, 15};

    // std::vector<int> ia = {0, 4, 7, 10, 14, 17};
    // std::vector<int> ja = {
    //     0, 1, 3, 4,
    //     0, 1, 3,
    //     2, 3, 4,
    //     0, 1, 2, 3,
    //     0, 2, 4};
    // std::vector<double> a = {
    //     11, 1, 3, 4,
    //     1, 12, 2,
    //     13, 1, 2,
    //     3, 2, 1, 14,
    //     4, 2, 15};
    //
    // std::vector<double> b = {23, 29, 17, 36, 21};
    // std::vector<double> diag = {11, 12, 13, 14, 15};

    std::vector<double> res(5);

    int nodes = 5;
    // std::vector<std::vector<double>> matrix(nodes, std::vector<double>(nodes, 0));
    // buildMatrix(ia, ja, a, matrix);
    // printMatrix(matrix);

    solve(ia, ja, a, b, diag, res);


    // for (int i = 0; i < 5; i++) {
    //     std::cout << res[i] << " ";
    // }
    //
    // std::cout << std::endl;


    std::vector<double> b_pred(5);
    std::vector<double> residual(5);
    spMV(ia, ja, a, res, b_pred);
    for (int i = 0; i < 5; i++) {
        residual[i] = b[i] - b_pred[i];
        std::cout << res[i] << std::endl;
    }

    double norm = dot(residual, residual);
    // std::cout << norm << std::endl;

    return norm > 1e-10;
}