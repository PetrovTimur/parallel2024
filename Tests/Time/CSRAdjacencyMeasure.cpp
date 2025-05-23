#include <iostream>
#include <tuple>
#include <vector>
#include "Utilities/input.h"
#include "Solver/csr.h"
#include "omp.h"
#include "Solver/Kernels/mathfunc.h"

int main() {
    const std::string elementsTxtPath = PROJECT_SOURCE_DIR "/Data/mesh100K/elements.txt";
    const std::string elementsDatPath = PROJECT_SOURCE_DIR "/Data/mesh100K/elements.dat";
    int Nn, Ne, M;
    int* rawData = nullptr;
    try {
        std::tie(Nn, Ne, M, rawData) = readData(elementsTxtPath, elementsDatPath);
    } catch (const std::runtime_error &e) {
        std::cerr << "Error reading data: " << e.what() << std::endl;
        return 1;
    }

    // Build element-to-node CSR (ia_en and ja_en) by filtering out -1 entries.
    std::vector<int> ia_en_vector;
    std::vector<int> ja_en_vector;
    ia_en_vector.push_back(0); // starting pointer
    for (int elem = 0; elem < Ne; elem++) {
        int count = 0;
        for (int j = 0; j < M; j++) {
            int node = rawData[elem * M + j];
            if (node != -1) {
                ja_en_vector.push_back(node);
                count++;
            }
        }
        ia_en_vector.push_back(ia_en_vector.back() + count);
    }
    // Free the raw data as it is no longer needed.
    delete[] rawData;

    // auto ia_en = new int[ia_en_vector.size()];
    // arrCopy(ia_en, ia_en_vector.data(), ia_en_vector.size());
    // auto ja_en = new int[ja_en_vector.size()];
    // arrCopy(ja_en, ja_en_vector.data(), ja_en_vector.size());

    // Transpose the CSR matrix to obtain node-to-element connectivity (ia_ne, ja_ne).
    // The transpose function treats the original CSR as a matrix with Ne rows and Nn columns.
    int *ia_ne, *ja_ne;
    double start_tr = omp_get_wtime();
    transposeCSR(ia_en_vector, ja_en_vector, Nn, ia_ne, ja_ne);
    double end_tr = omp_get_wtime();
    double duration_tr = (end_tr - start_tr) * 1000; // convert to milliseconds
    std::cout << "Time to transpose: " << duration_tr << " ms" << std::endl;


    // Measure time for building adjacency matrix using sets variant.
    double start_sets = omp_get_wtime();
    auto adj_sets = buildAdjacencyMatrixCSRUsingSets(ia_en_vector.data(), ja_en_vector.data(), ia_ne, ja_ne, Ne, Nn);
    double end_sets = omp_get_wtime();
    double duration_sets = (end_sets - start_sets) * 1000; // convert to milliseconds
    std::cout << "Time (using sets): " << duration_sets << " ms" << std::endl;

    // Measure time for building adjacency matrix using sort variant.
    double start_sort = omp_get_wtime();
    auto adj_sort = buildAdjacencyMatrixCSRUsingSort(ia_en_vector.data(), ja_en_vector.data(), ia_ne, ja_ne, Ne, Ne);
    double end_sort = omp_get_wtime();
    double duration_sort = (end_sort - start_sort) * 1000; // convert to milliseconds
    std::cout << "Time (using sort): " << duration_sort << " ms" << std::endl;

    // Measure time for building adjacency matrix using single CSR variant.
    double start_single = omp_get_wtime();
    auto adj_single = buildAdjacencyMatrixCSR(ia_ne, ja_ne, Ne, Nn);
    double end_single = omp_get_wtime();
    double duration_single = (end_single - start_single) * 1000; // convert to milliseconds
    std::cout << "Time (using single): " << duration_single << " ms" << std::endl;

    // Clean up dynamically allocated arrays.
//    delete[] ia_ne;
//    delete[] ja_ne;
    delete[] adj_sets.first;
    delete[] adj_sets.second;
    delete[] adj_sort.first;
    delete[] adj_sort.second;
    delete[] adj_single.first;
    delete[] adj_single.second;

    return 0;
}
