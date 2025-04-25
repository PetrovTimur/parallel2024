#include "csr.h"
#include "Kernels/mathfunc.h"

#include <algorithm>
#include <cmath>
#include <set>
#include <vector>
#include <utility>
#include <numeric>
#include <omp.h>
#include <unordered_map>
// #include <unordered_set>

#ifdef USE_MPI
void localizeCSR(const int *ia, const int size, int *ja, std::unordered_map<int, int> &G2L) {
    for (int i = 0; i < size - 1; i++) {
        for (int j = ia[i]; j < ia[i + 1]; j++) {
            ja[j] = G2L[ja[j]];
        }
    }
}

void constructG2L(std::vector<int> &ia_en, std::vector<int> &ja_en, std::unordered_map<int, int> &G2L) {
    int k = 0;
    for (int i = 0; i < ja_en.size(); i++) {
        if (!G2L.count(ja_en[i])) {
            G2L[ja_en[i]] = k++;
        }
    }
}

void fillG2L(std::vector<int> &L2G, std::unordered_map<int, int> &G2L) {
    for (int j = 0; j < L2G.size(); j++) {
        G2L[L2G[j]] = j;
    }
}

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
        diag[i] = 1. / a[k_i];
        // std::cout << "el: " << i << ", j: " << ja[k_i] << ", val = " << a[k_i] << std::endl;
    }

    #pragma omp parallel for proc_bind(master)
    for (int i = 0; i < b.size(); i++)
        b[i] = std::sin(L2G[i]);
}

void makeIncidenceMatrixCSR(const int Nx, int Ny, const int K1, const int K2, std::vector<int> &L2G, std::vector<int> &ia, std::vector<int> &ja, std::vector<int> &Part) {
    ia.push_back(0);

    std::vector<int> new_L2G;
    std::vector<int> new_Part;

    for (int i = 0; i < L2G.size(); i++) {
        int cell = L2G[i];
        int global_idx = cell + cell / (K1 + K2) * K2 + std::max(0, cell % (K1 + K2) - K1);

        if (cell % (K1 + K2) < K1) {
            ja.push_back(cell + cell / Nx);
            ja.push_back(cell + cell / Nx + 1);
            ja.push_back((cell + Nx + 1) + cell / Nx);
            ja.push_back((cell + Nx + 1) + cell / Nx + 1);
            ia.push_back(ia[ia.size() - 1] + 4);
            new_L2G.push_back(global_idx);
            new_Part.push_back(Part[i]);
        } else {
            ja.push_back(cell + cell / Nx);
            ja.push_back(cell + cell / Nx + 1);
            ja.push_back((cell + Nx + 1) + cell / Nx);
            ia.push_back(ia[ia.size() - 1] + 3);
            new_L2G.push_back(global_idx);
            new_Part.push_back(Part[i]);

            ja.push_back(cell + cell / Nx + 1);
            ja.push_back((cell + Nx + 1) + cell / Nx);
            ja.push_back((cell + Nx + 1) + cell / Nx + 1);
            ia.push_back(ia[ia.size() - 1] + 3);
            new_L2G.push_back(global_idx + 1);
            new_Part.push_back(Part[i]);
        }
    }

    L2G.clear();
    L2G.insert(L2G.end(), new_L2G.begin(), new_L2G.end());

    Part.clear();
    Part.insert(Part.end(), new_Part.begin(), new_Part.end());
}

void buildAdjacencyMatrixCSRUsingSort(const int *ia_en, const int *ja_en, const int *ia_ne, const int *ja_ne, std::vector<int> &ia_adj, std::vector<int> &
                                      ja_adj, const int Ne, std::vector<int> &Part, const int MyID) {
    ia_adj.push_back(0);

    for (int element1 = 0; element1 < Ne; ++element1) {
        if (Part[element1] != MyID)
            continue;

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
        ia_adj.push_back(ia_adj.back() + count);
    }
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


void makeIncidenceMatrixCSR(const int Nx, int Ny, const int K1, const int K2, const int i_start, const int i_end, const int j_start, const int j_end, std::vector<int> &ia, std::vector<int> &ja) {
    ia.push_back(0);

    for (int i = i_start; i <= i_end; i++) {
        for (int j = j_start; j <= j_end; j++) {
            const int cell = i * Nx + j;
            if (cell % (K1 + K2) < K1) {
                ja.push_back(cell + cell / Nx);
                ja.push_back(cell + cell / Nx + 1);
                ja.push_back((cell + Nx + 1) + cell / Nx);
                ja.push_back((cell + Nx + 1) + cell / Nx + 1);
                ia.push_back(ia[ia.size() - 1] + 4);
            } else {
                ja.push_back(cell + cell / Nx);
                ja.push_back(cell + cell / Nx + 1);
                ja.push_back((cell + Nx + 1) + cell / Nx);
                ia.push_back(ia[ia.size() - 1] + 3);

                ja.push_back(cell + cell / Nx + 1);
                ja.push_back((cell + Nx + 1) + cell / Nx);
                ja.push_back((cell + Nx + 1) + cell / Nx + 1);
                ia.push_back(ia[ia.size() - 1] + 3);
            }
        }
    }
}

template <typename T>
void fillCSR(const int *ia, const int *ja, T *a, T *b, T *diag, const int size) {
    #pragma omp parallel for default(none) shared(ia, ja, a, b, diag, size)
    for (int i = 0; i < size; i++) {
        T sum = 0;
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
        diag[i] = 1. / a[k_i];
        // std::cout << "el: " << i << ", j: " << ja[k_i] << ", val = " << a[k_i] << std::endl;
    }

    #pragma omp parallel for default(none) shared(b, size)
    for (int i = 0; i < size; i++)
        b[i] = std::sin(i);
}

#ifdef USE_CUDA
template void fillCSR<float>(const int *ia, const int *ja, float *a, float *b, float *diag, int size);
#else
template void fillCSR<double>(const int *ia, const int *ja, double *a, double *b, double *diag, int size);
#endif

std::pair<int*, int*> buildAdjacencyMatrixCSRUsingSort(const int *ia_en, const int *ja_en, const int *ia_ne, const int *ja_ne, const int Ne, const int Nn) {
    int num_threads;
    #pragma omp parallel default(none) shared(num_threads)
    {
        #pragma omp single
        {
            num_threads = omp_get_num_threads();
        }
    }

    auto ia_adj = new int[Ne + 1];
    ia_adj[0] = 0;
    std::vector<std::vector<int>> ja_local(num_threads);
    std::vector<int> ia_local(num_threads + 1);

    #pragma omp parallel default(none) shared(ia_en, ja_en, ia_ne, ja_ne, ia_adj, ia_local, ja_local, num_threads, Ne)
    {
        int tid = omp_get_thread_num();
        int chunk = (Ne + num_threads - 1) / num_threads;
        int start = tid * chunk;
        int end   = std::min(start + chunk, Ne);
        for (int element1 = start; element1 < end; ++element1) {
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

            ja_local[tid].push_back(adjacent[0]);
            int count = 1;
            for (unsigned int i = 1; i < adjacent.size(); ++i) {
                if (adjacent[i] != adjacent[i - 1]) {
                    ja_local[tid].push_back(adjacent[i]);
                    count++;
                }
            }
            ia_adj[element1 + 1] = count;
            ia_local[tid + 1] += count;
        }
    }

    scan(ia_adj, ia_adj, Ne + 1);

    for (int tid = 1; tid <= num_threads; ++tid) {
        ia_local[tid] += ia_local[tid - 1];
    }

    auto ja_adj = new int[ia_local[num_threads]];

    #pragma omp parallel for default(none) shared(ja_adj, ia_local, ja_local, num_threads)
    for (int tid = 0; tid < num_threads; tid++) {
        for (int i = ia_local[tid]; i < ia_local[tid + 1]; i++) {
            ja_adj[i] = ja_local[tid][i - ia_local[tid]];
        }
    }

    return std::make_pair(ia_adj, ja_adj);
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
#endif

void buildMatrixFromCSR(std::vector<int> &ia, std::vector<int> &ja, std::vector<std::vector<bool>> &matrix) {
    for (unsigned int i = 0; i < ia.size() - 1; i++) {
        for (int k = ia[i]; k < ia[i + 1]; k++) {
            matrix[i][ja[k]] = true;
        }
    }
}

void buildFilledMatrix(std::vector<int> &ia, std::vector<int> &ja, std::vector<double> &a,
    std::vector<std::vector<double>> &matrix) {
    for (unsigned int i = 0; i < ia.size() - 1; i++) {
        for (int k = ia[i]; k < ia[i + 1]; k++) {
            matrix[i][ja[k]] = a[k];
        }
    }
}

void transposeCSR(std::vector<int> &ia, std::vector<int> &ja, const int Nn, int *&ia_new, int *&ja_new) {
    const int nnz = ia[ia.size() - 1];
    ia_new = new int[Nn + 1];
    ja_new = new int[nnz];

    arrInit(ia_new, 0, Nn + 1);

    #pragma omp parallel for default(none) shared(ia_new, ja, nnz)
    for (int i = 0; i < nnz; ++i) {
        #pragma omp atomic update
        ia_new[ja[i] + 1]++;
    }

    scan(ia_new, ia_new, Nn + 1);

    const auto buf = new int[Nn + 1];
    arrCopy(buf, ia_new, Nn + 1);

    #pragma omp parallel for default(none) shared(ia, ja, ja_new, buf)
    for (unsigned int i = 0; i < ia.size() - 1; ++i) {
        for (int k = ia[i]; k < ia[i + 1]; ++k) {
            const int col = ja[k];
            int dest_pos;
            #pragma omp atomic capture
            dest_pos = buf[col]++;
            ja_new[dest_pos] = i;
        }
    }

    delete[] buf;
}