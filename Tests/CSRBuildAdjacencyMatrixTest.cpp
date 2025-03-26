#include "Solver/csr.h"
#include <iostream>
#include <vector>

// Test adjacency matrix from incidence matrix with multiple connections
bool IncidenceToAdjacencyMultipleConnections() {
    // Incidence matrix: Elements as rows, nodes as columns
    int ia[] = {0, 2, 3, 5};  // Row pointers for elements
    int ja[] = {0, 1, 2, 0, 2};  // Column indices (nodes)
    int rows = 3;  // Number of elements
    int cols = 3;  // Number of nodes

    // First transpose to get nodes as rows
    auto transposed = transposeCSR(ia, ja, rows, cols);

    // Now build adjacency matrix
    auto result = buildAdjacencyMatrixCSR(transposed.first, transposed.second, rows, cols);

    // Expected adjacency matrix (node->node connections)
    int expected_ia[] = {0, 2, 4, 7};
    int expected_ja[] = {0, 2, 1, 2, 0, 1, 2};

    bool success = true;

    // Check ia values
    for (int i = 0; i <= rows; i++) {
        if (result.first[i] != expected_ia[i]) {
            std::cout << "Error in ia[" << i << "]: expected " << expected_ia[i]
                      << ", got " << result.first[i] << std::endl;
            success = false;
        }
    }

    // Check ja values
    for (int i = 0; i < result.first[rows]; i++) {
        if (result.second[i] != expected_ja[i]) {
            std::cout << "Error in ja[" << i << "]: expected " << expected_ja[i]
                      << ", got " << result.second[i] << std::endl;
            success = false;
        }
    }

    delete[] result.first;
    delete[] result.second;
    delete[] transposed.first;
    delete[] transposed.second;

    return success;
}

// Test simple path graph incidence to adjacency conversion
bool PathGraphIncidenceToAdjacency() {
    // Incidence matrix for path graph with 3 nodes and 2 edges
    int ia[] = {0, 2, 4};  // Row pointers for elements
    int ja[] = {0, 1, 1, 2};  // Column indices (nodes)
    int rows = 2;  // Number of elements
    int cols = 3;  // Number of nodes

    // First transpose to get nodes as rows
    auto transposed = transposeCSR(ia, ja, rows, cols);

    // Now build adjacency matrix
    auto result = buildAdjacencyMatrixCSR(transposed.first, transposed.second, rows, cols);

    // Expected adjacency matrix
    int expected_ia[] = {0, 2, 4};
    int expected_ja[] = {0, 1, 0, 1};

    bool success = true;

    // Check ia values
    for (int i = 0; i <= rows; i++) {
        if (result.first[i] != expected_ia[i]) {
            std::cout << "Error in ia[" << i << "]: expected " << expected_ia[i]
                      << ", got " << result.first[i] << std::endl;
            success = false;
        }
    }

    // Check ja values
    for (int i = 0; i < result.first[rows]; i++) {
        if (result.second[i] != expected_ja[i]) {
            std::cout << "Error in ja[" << i << "]: expected " << expected_ja[i]
                      << ", got " << result.second[i] << std::endl;
            success = false;
        }
    }

    delete[] result.first;
    delete[] result.second;
    delete[] transposed.first;
    delete[] transposed.second;

    return success;
}

// Test isolated nodes in incidence matrix
bool IsolatedNodesIncidenceToAdjacency() {
    // Incidence matrix with isolated node
    int ia[] = {0, 2};  // Row pointers for elements
    int ja[] = {0, 2};  // Column indices (nodes) - node 1 is isolated
    int rows = 1;  // Number of elements
    int cols = 3;  // Number of nodes

    // First transpose to get nodes as rows
    auto transposed = transposeCSR(ia, ja, rows, cols);

    // Now build adjacency matrix
    auto result = buildAdjacencyMatrixCSR(transposed.first, transposed.second, rows, cols);

    // Expected adjacency matrix (isolated node has no connections)
    int expected_ia[] = {0, 1};
    int expected_ja[] = {0};

    bool success = true;

    // Check ia values
    for (int i = 0; i <= rows; i++) {
        if (result.first[i] != expected_ia[i]) {
            std::cout << "Error in ia[" << i << "]: expected " << expected_ia[i]
                      << ", got " << result.first[i] << std::endl;
            success = false;
        }
    }

    // Check ja values
    for (int i = 0; i < result.first[rows]; i++) {
        if (result.second[i] != expected_ja[i]) {
            std::cout << "Error in ja[" << i << "]: expected " << expected_ja[i]
                      << ", got " << result.second[i] << std::endl;
            success = false;
        }
    }

    delete[] result.first;
    delete[] result.second;
    delete[] transposed.first;
    delete[] transposed.second;

    return success;
}

// Test completely disconnected graph
bool DisconnectedGraphIncidenceToAdjacency() {
    // Incidence matrix with no connections
    int ia[] = {0, 0, 0};  // Row pointers for elements - no elements
    int ja[] = {};  // No connections
    int rows = 2;  // Number of elements (dummy)
    int cols = 3;  // Number of nodes

    // First transpose to get nodes as rows
    auto transposed = transposeCSR(ia, ja, rows, cols);

    // Now build adjacency matrix
    auto result = buildAdjacencyMatrixCSR(transposed.first, transposed.second, rows, cols);

    // Expected adjacency matrix (all zeros)
    int expected_ia[] = {0, 0, 0};
    int expected_ja[] = {};

    bool success = true;

    // Check ia values
    for (int i = 0; i <= rows; i++) {
        if (result.first[i] != expected_ia[i]) {
            std::cout << "Error in ia[" << i << "]: expected " << expected_ia[i]
                      << ", got " << result.first[i] << std::endl;
            success = false;
        }
    }

    delete[] result.first;
    delete[] result.second;
    delete[] transposed.first;
    delete[] transposed.second;

    return success;
}

int main() {
    bool all_passed = true;

    std::cout << "Testing incidence to adjacency with multiple connections... ";
    all_passed &= IncidenceToAdjacencyMultipleConnections();
    std::cout << (all_passed ? "PASSED\n" : "FAILED\n");

    std::cout << "Testing path graph incidence to adjacency... ";
    all_passed &= PathGraphIncidenceToAdjacency();
    std::cout << (all_passed ? "PASSED\n" : "FAILED\n");

    std::cout << "Testing incidence to adjacency with isolated nodes... ";
    all_passed &= IsolatedNodesIncidenceToAdjacency();
    std::cout << (all_passed ? "PASSED\n" : "FAILED\n");

    std::cout << "Testing disconnected graph incidence to adjacency... ";
    all_passed &= DisconnectedGraphIncidenceToAdjacency();
    std::cout << (all_passed ? "PASSED\n" : "FAILED\n");

    return all_passed ? 0 : 1;
}