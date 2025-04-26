#include "Solver/csr.h"
#include <iostream>
#include <vector>

// Test transpose of a square matrix
bool TestSquareMatrixTranspose() {
    std::vector<int> ia = {0, 2, 3, 5};  // Row pointers
    std::vector<int> ja = {0, 1, 2, 0, 2};  // Column indices
    int rows = 3;
    int cols = 3;

    int *ia_ne, *ja_ne;
    transposeCSR(ia, ja, cols, ia_ne, ja_ne);
    auto result = std::make_pair(ia_ne, ja_ne);

    int expected_ia[] = {0, 2, 3, 5};
    int expected_ja[] = {0, 2, 0, 1, 2};

    bool success = true;

    // Check ia values
    for (int i = 0; i <= cols; i++) {
        if (result.first[i] != expected_ia[i]) {
            std::cout << "Error in ia[" << i << "]: expected " << expected_ia[i]
                      << ", got " << result.first[i] << std::endl;
            success = false;
        }
    }

    // Check ja values
    for (int i = 0; i < ia[rows]; i++) {
        if (result.second[i] != expected_ja[i]) {
            std::cout << "Error in ja[" << i << "]: expected " << expected_ja[i]
                      << ", got " << result.second[i] << std::endl;
            success = false;
        }
    }

    delete[] result.first;
    delete[] result.second;

    return success;
}

// Test transpose of a rectangular matrix
bool TestRectangularMatrixTranspose() {
    std::vector<int> ia = {0, 2, 4};  // Row pointers
    std::vector<int> ja = {0, 2, 1, 2};  // Column indices
    int rows = 2;
    int cols = 3;

    int *ia_ne, *ja_ne;
    transposeCSR(ia, ja, cols, ia_ne, ja_ne);
    auto result = std::make_pair(ia_ne, ja_ne);

    int expected_ia[] = {0, 1, 2, 4};
    int expected_ja[] = {0, 1, 0, 1};

    bool success = true;

    // Check ia values
    for (int i = 0; i <= cols; i++) {
        if (result.first[i] != expected_ia[i]) {
            std::cout << "Error in ia[" << i << "]: expected " << expected_ia[i]
                      << ", got " << result.first[i] << std::endl;
            success = false;
        }
    }

    // Check ja values
    for (int i = 0; i < ia[rows]; i++) {
        if (result.second[i] != expected_ja[i]) {
            std::cout << "Error in ja[" << i << "]: expected " << expected_ja[i]
                      << ", got " << result.second[i] << std::endl;
            success = false;
        }
    }

    delete[] result.first;
    delete[] result.second;

    return success;
}

// Test transpose of an empty matrix
bool TestEmptyMatrixTranspose() {
    std::vector<int> ia = {0, 0};  // Empty matrix
    std::vector<int> ja = {};      // No elements
    // int rows = 1;
    int cols = 1;

    int *ia_ne, *ja_ne;
    transposeCSR(ia, ja, cols, ia_ne, ja_ne);
    auto result = std::make_pair(ia_ne, ja_ne);

    bool success = true;

    // Transposed matrix should also be empty
    if (result.first[0] != 0 || result.first[1] != 0) {
        std::cout << "Error in empty matrix transpose: expected ia[0]=0, ia[1]=0, got "
                  << result.first[0] << ", " << result.first[1] << std::endl;
        success = false;
    }

    delete[] result.first;
    delete[] result.second;

    return success;
}

// Test transpose of a diagonal matrix
bool TestDiagonalMatrixTranspose() {
    std::vector<int> ia = {0, 1, 2, 3};
    std::vector<int> ja = {0, 1, 2};
    int rows = 3;
    int cols = 3;

    int *ia_ne, *ja_ne;
    transposeCSR(ia, ja, cols, ia_ne, ja_ne);
    auto result = std::make_pair(ia_ne, ja_ne);

    int expected_ia[] = {0, 1, 2, 3};
    int expected_ja[] = {0, 1, 2};

    bool success = true;

    // Check ia values
    for (int i = 0; i <= cols; i++) {
        if (result.first[i] != expected_ia[i]) {
            std::cout << "Error in ia[" << i << "]: expected " << expected_ia[i]
                      << ", got " << result.first[i] << std::endl;
            success = false;
        }
    }

    // Check ja values
    for (int i = 0; i < ia[rows]; i++) {
        if (result.second[i] != expected_ja[i]) {
            std::cout << "Error in ja[" << i << "]: expected " << expected_ja[i]
                      << ", got " << result.second[i] << std::endl;
            success = false;
        }
    }

    delete[] result.first;
    delete[] result.second;

    return success;
}

int main() {
    bool all_passed = true;

    std::cout << "Testing square matrix transpose... ";
    all_passed &= TestSquareMatrixTranspose();
    std::cout << (all_passed ? "PASSED\n" : "FAILED\n");

    std::cout << "Testing rectangular matrix transpose... ";
    all_passed &= TestRectangularMatrixTranspose();
    std::cout << (all_passed ? "PASSED\n" : "FAILED\n");

    std::cout << "Testing empty matrix transpose... ";
    all_passed &= TestEmptyMatrixTranspose();
    std::cout << (all_passed ? "PASSED\n" : "FAILED\n");

    std::cout << "Testing diagonal matrix transpose... ";
    all_passed &= TestDiagonalMatrixTranspose();
    std::cout << (all_passed ? "PASSED\n" : "FAILED\n");

    return all_passed ? 0 : 1;
}