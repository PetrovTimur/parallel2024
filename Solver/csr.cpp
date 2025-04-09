#include "csr.h"
#include "Kernels/mathfunc.h"

#include <algorithm>
#include <cmath>
#include <set>
#include <vector>
#include <utility>
#include <numeric>
// #include <unordered_set>

#ifdef USE_MPI
void makeCSR(int Nx, int Ny, int K1, int K2, int i_start, int i_end, int j_start, int j_end, std::vector<int> &G2L, std::vector<int> &ia, std::vector<int> &ja) {
    int cycle_length = K1 + K2;
    // std::vector<int> extras(i_end - i_start + 1);

    // #pragma omp parallel for proc_bind(spread)
    // for (int i = i_start; i <= i_end; i++) {
    //     extras[i - i_start] = i * Nx / cycle_length * K2 + std::max(0, i * Nx % cycle_length - K1);
        // std::cout << extras[i] << " ";
    // }
    // std::cout << std::endl;

    std::vector<std::vector<int>> temp_ia(i_end - i_start + 1, std::vector<int>(j_end - j_start + 1));
    std::vector<std::vector<int>> temp_ja(i_end - i_start + 1, std::vector<int>(7 * (j_end - j_start + 1)));

    // int k = 0;
    #pragma omp parallel for proc_bind(master)
    for (int i = i_start; i <= i_end; i++) {
        // std::cout << i << " " << k << std::endl;
        // int k =
        // int k = 2 * Nx * i + (Nx + 1) * std::max(0, 2 * i - 1) + i * (Nx + 1) + (i > 0 ? extras[i] + extras[i - 1] : 0);
        // std::cout << i << " " << k << std::endl;
        int temp_k = 0;
        for (int j = j_start; j <= j_end; j++) {
            int ttemp_k = temp_k;;
            int flat_idx = i * (Nx + 1) + j;
            // std::cout << i << " " << j << " " << k << " " << flat_idx << " " << G2L[flat_idx] << std::endl;
            // int top_idx = (i - 1) * (Nx + 1) + j;
            // int right_idx = i * (Nx + 1) + j + 1;
            // int bottom_idx = (i + 1) * (Nx + 1) + j;
            // int left_idx = i * (Nx + 1) + j - 1;
            // ia[G2L[flat_idx]] = k;

            if (i != 0) {
                // ja[k++] = G2L[(i - 1) * (Nx + 1) + j];
                temp_ja[i - i_start][temp_k++] = G2L[(i - 1) * (Nx + 1) + j];
            }
            if (i != 0 && j != Nx && ((i - 1) * Nx + j) % cycle_length >= K1) {
                // ja[k++] = G2L[(i - 1) * (Nx + 1) + j + 1];
                temp_ja[i - i_start][temp_k++] = G2L[(i - 1) * (Nx + 1) + j + 1];
            }
            if (j != 0) {
                // ja[k++] = G2L[i * (Nx + 1) + j - 1];
                temp_ja[i - i_start][temp_k++] = G2L[i * (Nx + 1) + j - 1];
            }

            // ja[k++] = G2L[flat_idx];
            temp_ja[i - i_start][temp_k++] = G2L[flat_idx];

            if (j != Nx) {
                // ja[k++] = G2L[i * (Nx + 1) + j + 1];
                temp_ja[i - i_start][temp_k++] = G2L[i * (Nx + 1) + j + 1];
            }
            if  (j != 0 && i != Ny && (i * Nx + j - 1) % cycle_length >= K1) {
                // ja[k++] = G2L[(i + 1) * (Nx + 1) + j - 1];
                temp_ja[i - i_start][temp_k++] = G2L[(i + 1) * (Nx + 1) + j - 1];
            }
            if (i != Ny) {
                // ja[k++] = G2L[(i + 1) * (Nx + 1) + j];
                temp_ja[i - i_start][temp_k++] = G2L[(i + 1) * (Nx + 1) + j];
            }

            temp_ia[i - i_start][j - j_start] = temp_k - ttemp_k;
        }
        temp_ja[i - i_start].resize(temp_k);
    }
    // printMatrix(temp_ia);
    // std::vector<int> test_ia = {0};
    // std::vector<int> test_ja;

    for (const auto& row : temp_ia) {
        ia.insert(ia.end(), row.begin(), row.end());
    }

    for (const auto& row : temp_ja) {
        ja.insert(ja.end(), row.begin(), row.end());
    }


    // ia[ia.size() - 1] = k;
    // ja.resize(k);

    std::partial_sum(ia.begin(), ia.end(), ia.begin());
    // std::cout << "Test IA: " << test_ia.size() << std::endl;
    // printVector(test_ia);
    // std::cout << "IA: " << ia.size() << std::endl;
    // printVector(ia);
    // std::cout << "Test JA: " << test_ja.size() << std::endl;
    // printVector(test_ja);
    // std::cout << "JA: " << ja.size() << std::endl;
    // printVector(ja);
}

#else
void makeCSR(int Nx, int Ny, int K1, int K2, std::vector<int> &ia, std::vector<int> &ja) {
    int cycle_length = K1 + K2;
    std::vector<int> extras(Ny + 1);
    #pragma omp parallel for proc_bind(spread)
    for (int i = 0; i < Ny + 1; i++) {
        extras[i] = i * Nx / cycle_length * K2 + std::max(0, i * Nx % cycle_length - K1);
        // std::cout << extras[i] << " ";
    }
    // std::cout << std::endl;

    // int k = 0;
    #pragma omp parallel for proc_bind(spread)
    for (int i = 0; i < Ny + 1; i++) {
        // std::cout << i << " " << k << std::endl;
        int k = 2 * Nx * i + (Nx + 1) * std::max(0, 2 * i - 1) + i * (Nx + 1) + (i > 0 ? extras[i] + extras[i - 1] : 0);
        // std::cout << i << " " << k << std::endl;
        for (int j = 0; j < Nx + 1; j++) {
            int flat_idx = i * (Nx + 1) + j;
            // int top_idx = (i - 1) * (Nx + 1) + j;
            // int right_idx = i * (Nx + 1) + j + 1;
            // int bottom_idx = (i + 1) * (Nx + 1) + j;
            // int left_idx = i * (Nx + 1) + j - 1;
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
    ia[(Nx + 1) * (Ny + 1)] = ja.size();
}

#endif
std::pair<int *, int *> makeIncidenceMatrixCSR(int Nx, int Ny, int K1, int K2, int Ne, int nnz) {
    auto ia = new int[Ne + 1];
    auto ja = new int[nnz];
    ia[0] = 0;

    int offset = 0;
    for (int cell = 0; cell < Nx * Ny; cell++) {
        if (cell % (K1 + K2) < K1) {
            int element_number = cell + offset;
            ia[element_number + 1] = ia[element_number] + 4;

            ja[ia[element_number]] = cell + cell / Nx;
            ja[ia[element_number] + 1] = cell + cell / Nx + 1;
            ja[ia[element_number] + 2] = (cell + Nx + 1) + cell / Nx;
            ja[ia[element_number] + 3] = (cell + Nx + 1) + cell / Nx + 1;
        } else {
            int element_number = cell + offset;
            ia[element_number + 1] = ia[element_number] + 3;
            ja[ia[element_number]] = cell + cell / Nx;
            ja[ia[element_number] + 1] = cell + cell / Nx + 1;
            ja[ia[element_number] + 2] = (cell + Nx + 1) + cell / Nx;
            // ja[ia[element_number] + 3] = (cell + Nx + 1) + cell / Nx + 1;

            element_number = cell + ++offset;
            ia[element_number + 1] = ia[element_number] + 3;
            // ja[ia[element_number]] = cell + cell / Nx;
            ja[ia[element_number]] = cell + cell / Nx + 1;
            ja[ia[element_number] + 1] = (cell + Nx + 1) + cell / Nx;
            ja[ia[element_number] + 2] = (cell + Nx + 1) + cell / Nx + 1;

            // offset++;
        }
    }

    return std::make_pair(ia, ja);
}

#ifdef USE_MPI
void fillCSR(std::vector<int> &ia, std::vector<int> &ja, std::vector<int> &L2G, std::vector<double> &a, std::vector<double> &b,
    std::vector<double> &diag) {
    #pragma omp parallel for proc_bind(master)
    for (int i = 0; i < ia.size() - 1; i++) {
        double sum = 0;
        int k_i = 0;
        for (int k = ia[i]; k < ia[i + 1]; k++) {
            int j = ja[k];
            if (i != j) {
                a[k] = std::cos(L2G[i] * L2G[j] + L2G[i] + L2G[j]);
                sum += std::abs(a[k]);
            } else
                k_i = k;
            // std::cout << "el: " << i << ", j: " << j << ", val = " << a[k] << std::endl;
        }

        a[k_i] = 1.234 * sum;
        diag[i] = a[k_i];
        // std::cout << "el: " << i << ", j: " << ja[k_i] << ", val = " << a[k_i] << std::endl;
    }

    #pragma omp parallel for proc_bind(master)
    for (int i = 0; i < b.size(); i++)
        b[i] = std::sin(L2G[i]);
}
#else
void fillCSR(const int *ia, const int *ja, double *a, double *b, double *diag, const int size) {
    #pragma omp parallel for proc_bind(spread)
    for (int i = 0; i < size; i++) {
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

    #pragma omp parallel for proc_bind(spread)
    for (int i = 0; i < size; i++)
        b[i] = std::sin(i);
}
#endif

void buildMatrixFromCSR(std::vector<int> &ia, std::vector<int> &ja, std::vector<std::vector<bool>> &matrix) {
    for (int i = 0; i < ia.size() - 1; i++) {
        for (int k = ia[i]; k < ia[i + 1]; k++) {
            matrix[i][ja[k]] = true;
        }
    }
}

void buildFilledMatrix(std::vector<int> &ia, std::vector<int> &ja, std::vector<double> &a,
    std::vector<std::vector<double>> &matrix) {
    for (int i = 0; i < ia.size() - 1; i++) {
        for (int k = ia[i]; k < ia[i + 1]; k++) {
            matrix[i][ja[k]] = a[k];
        }
    }
}

std::pair<int *, int *> transposeCSR(const int *ia, const int *ja, const int Ne, const int Nn) {
    const int nnz = ia[Ne];
    auto ia_new = new int[Nn + 1];
    auto ja_new = new int[nnz];

    // Zero out ia_new
    for (int i = 0; i <= Nn; i++) {
        ia_new[i] = 0;
    }

    // Count the number of entries in each column
    for (int i = 0; i < nnz; ++i) {
        ia_new[ja[i] + 1]++;
    }

    // Compute the cumulative sum to get the new row pointers
    for (int i = 1; i <= Nn; ++i) {
        ia_new[i] += ia_new[i - 1];
    }

    const auto buf = new int[Nn + 1];

    // Zero out buf
    for (int i = 0; i <= Nn; ++i) {
        buf[i] = ia_new[i];
    }

    // Fill the column indices
    for (int i = 0; i < Ne; ++i) {
        for (int k = ia[i]; k < ia[i + 1]; ++k) {
            const int col = ja[k];
            const int dest_pos = buf[col]++;
            ja_new[dest_pos] = i;
        }
    }

    delete[] buf;

    return std::make_pair(ia_new, ja_new);
}

std::pair<int *, int *> buildAdjacencyMatrixCSR(const int *ia_ne, const int *ja_ne, const int Ne, const int Nn) {
    auto adjacency_list = new std::set<int>[Ne];

    // Build adjacency list for elements
    for (int node = 0; node < Nn; ++node) {
        for (int k = ia_ne[node]; k < ia_ne[node + 1]; ++k) {
            const int element1 = ja_ne[k];
            for (int l = ia_ne[node]; l < ia_ne[node + 1]; ++l) {
                const int element2 = ja_ne[l];
                adjacency_list[element1].insert(element2);
            }
        }
    }

    // Convert adjacency list to CSR format
    int nnz = 0;
    for (int i = 0; i < Ne; ++i) {
        nnz += adjacency_list[i].size();
    }

    auto ia_adj = new int[Ne + 1];
    auto ja_adj = new int[nnz];

    ia_adj[0] = 0;
    int index = 0;
    for (int i = 0; i < Ne; ++i) {
        ia_adj[i + 1] = ia_adj[i] + adjacency_list[i].size();
        std::vector<int> sorted_neighbors(adjacency_list[i].begin(), adjacency_list[i].end());
        // std::sort(sorted_neighbors.begin(), sorted_neighbors.end());
        for (const int neighbor : sorted_neighbors) {
            ja_adj[index++] = neighbor;
        }
    }

    delete[] adjacency_list;

    return std::make_pair(ia_adj, ja_adj);
}

std::pair<int *, int *> buildAdjacencyMatrixCSRUsingSets(const int *ia_en, const int *ja_en, const int *ia_ne, const int *ja_ne, const int Ne, const int Nn) {
    std::vector<int> ia_adj(Ne + 1);
    std::vector<int> ja_adj;

    for (int element1 = 0; element1 < Ne; ++element1) {
        std::set<int> adjacent;
        for (int k = ia_en[element1]; k < ia_en[element1 + 1]; ++k) {
            const int node = ja_en[k];

            for (int l = ia_ne[node]; l < ia_ne[node + 1]; ++l) {
                const int element2 = ja_ne[l];
                adjacent.insert(element2);
            }
        }
        adjacent.insert(element1);

        ia_adj[element1 + 1] = ia_adj[element1] + adjacent.size();
        for (auto element: adjacent) {
            ja_adj.push_back(element);
        }
    }

    auto ia_adj_new = new int[ia_adj.size()];
    auto ja_adj_new = new int[ja_adj.size()];
    arrCopy(ia_adj_new, ia_adj.data(), ia_adj.size());
    arrCopy(ja_adj_new, ja_adj.data(), ja_adj.size());

    return std::make_pair(ia_adj_new, ja_adj_new);
}


std::pair<int*, int*> buildAdjacencyMatrixCSRUsingSort(const int *ia_en, const int *ja_en, const int *ia_ne, const int *ja_ne, const int Ne, const int Nn) {
    std::vector<int> ia_adj(Ne + 1);
    std::vector<int> ja_adj;

    for (int element1 = 0; element1 < Ne; ++element1) {
        std::vector<int> adjacent;
        for (int k = ia_en[element1]; k < ia_en[element1 + 1]; ++k) {
            const int node = ja_en[k];

            for (int l = ia_ne[node]; l < ia_ne[node + 1]; ++l) {
                const int element2 = ja_ne[l];
                adjacent.push_back(element2);
            }
        }
        adjacent.push_back(element1);

        std::sort(adjacent.begin(), adjacent.end());

        ja_adj.push_back(adjacent[0]);
        int count = 1;
        for (int i = 1; i < adjacent.size(); ++i) {
            if (adjacent[i] != adjacent[i - 1]) {
                ja_adj.push_back(adjacent[i]);
                count++;
            }
        }
        ia_adj[element1 + 1] = ia_adj[element1] + count;
    }

    auto ia_adj_new = new int[ia_adj.size()];
    auto ja_adj_new = new int[ja_adj.size()];
    arrCopy(ia_adj_new, ia_adj.data(), ia_adj.size());
    arrCopy(ja_adj_new, ja_adj.data(), ja_adj.size());

    return std::make_pair(ia_adj_new, ja_adj_new);
}