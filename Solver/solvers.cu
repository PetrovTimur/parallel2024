#include "cuda_runtime.h"
#include "thrust/device_vector.h"
#include "Solver/Kernels/mathfunc.cuh"
#include <Utilities/logger.h>
#include "solvers.cuh"
#include "omp.h"


#define checkCudaErrors(val) check((val), #val, __FILE__, __LINE__)
void check(cudaError_t err, const char* const func, const char* const file,
           const int line)
{
    if (err != cudaSuccess)
    {
        std::cerr << "CUDA Runtime Error at: " << file << ":" << line
                  << std::endl;
        std::cerr << cudaGetErrorString(err) << " " << func << std::endl;
        // We don't exit when we encounter CUDA errors in this example.
        std::exit(EXIT_FAILURE);
    }
}

#define getLastCudaError(msg) checkLast(msg, __FILE__, __LINE__)

inline void checkLast(const char *errorMessage, const char *file,
                               const int line) {
    cudaError_t err = cudaGetLastError();

    if (cudaSuccess != err) {
        fprintf(stderr,
                "%s(%i) : getLastCudaError() CUDA error :"
                " %s : (%d) %s.\n",
                file, line, errorMessage, static_cast<int>(err),
                cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
}

int solve(const int *ia, const int *ja, const float *a, const float *b, const float *diag, int size, float *res, const double eps, const int maxit) {
    const int N = size - 1;
    const int blocks = 128;
    const int threads = 256;

    thrust::device_vector<float> z(N);
    thrust::device_vector<float> p(N);
    thrust::device_vector<float> q(N);
    thrust::device_vector<float> x(N);
    thrust::device_vector<float> r(b, b + N);
    thrust::device_vector<float> diag_gpu(diag, diag + N);

    thrust::device_vector<float> rho(2);
    thrust::device_vector<int> ia_gpu(ia, ia + N + 1);
    thrust::device_vector<int> ja_gpu(ja, ja + ia[N]);
    thrust::device_vector<float> a_gpu(a, a + ia[N]);

    thrust::device_vector<float> vec_buf(blocks);

    float *buf_gpu, *norm_gpu;
    cudaMalloc(&buf_gpu, sizeof(float));
    cudaMalloc(&norm_gpu, sizeof(float));
    int k = 0;
    float buf, norm;

    // while (rho[1] > eps * eps && k < maxit) {
    //     if (k > 1) {
    //         beta = rho[1] / rho[0];
    //     } else {
    //         p = z;cuda-memcheck ./bin/CGSolver --Nx=4 --Ny=4 --K1=2 --K2=3
    //     }
    //
    //
    //
    //     k++;
    // }

    std::cout << "Starting: ";

    do {
        k++;

        LOG_INFO << "Iteration " << k << std::endl;
        std::cout << "Iteration " << k << std::endl;

        const auto start = omp_get_wtime();

        multiply<<<blocks, threads>>>(r.data().get(), diag_gpu.data().get(), z.data().get(), N);
        getLastCudaError("Kernel execution failed");

        rho[0] = rho[1];

        // std::cout << r.size() << std::endl << z.size() << std::endl << N << std::endl;

        assert(r.size() == z.size() && r.size() == N && "Size mismatch");
        dot_gpu(threads, blocks, r.data().get(), z.data().get(), vec_buf.data().get(), rho.data().get() + 1, N);
        getLastCudaError("Kernel execution failed");

        // rho[1] = *buf_gpu;

        if (k == 1)
            p = z;
        else {
            const float beta = rho[1] / rho[0];
            // TODO: move to GPU (?)
            // axpy(beta, p, z, N, p);
            axpy<<<blocks, threads>>>(beta, p.data().get(), z.data().get(), N, p.data().get());
            getLastCudaError("Kernel execution failed");

        }

        spMV<<<blocks, threads>>>(ia_gpu.data().get(), ja_gpu.data().get(), a_gpu.data().get(), p.data().get(), q.data().get(), N);
        getLastCudaError("Kernel execution failed");

        assert(p.size() == q.size() && p.size() == N && "Size mismatch");
        dot_gpu(threads, blocks, p.data().get(), q.data().get(), vec_buf.data().get(), buf_gpu, N);
        getLastCudaError("Kernel execution failed");

        cudaMemcpy(&buf, buf_gpu, sizeof(float), cudaMemcpyDeviceToHost);
        float alpha = rho[1] / buf;

        axpy<<<blocks, threads>>>(alpha, p.data().get(), x.data().get(), N, x.data().get());
        getLastCudaError("Kernel execution failed");

        axpy<<<blocks, threads>>>(-alpha, q.data().get(), r.data().get(), N, r.data().get());
        getLastCudaError("Kernel execution failed");


        LOG_INFO << "Time " << omp_get_wtime() - start << std::endl;
        // LOG_INFO << "rho = " << rho[0] << ", " << rho[1] << ", alpha = " << alpha << std::endl;
        dot_gpu(threads, blocks, x.data().get(), x.data().get(), vec_buf.data().get(), norm_gpu, N);
        // dot<<<blocks, threads>>>(x.data().get(), x.data().get(), norm_gpu, N);
        cudaMemcpy(&norm, norm_gpu, sizeof(float), cudaMemcpyDeviceToHost);
        LOG_INFO << "L2 norm: " << std::sqrt(norm) << std::endl;
        LOG << "--------------------------------------------------" << std::endl << std::endl;

    }
    while (rho[1] > eps * eps && k < maxit);

    // #pragma omp parallel for default(none) shared(res, x, N)
    cudaMemcpy(res, x.data().get(), N * sizeof(float), cudaMemcpyDeviceToHost);
    std::cout << x[0] << x[1] << std::endl;

    // delete[] z;
    // delete[] p;
    // delete[] q;
    // delete[] x;
    // delete[] r;
    // delete[] rho;

    cudaFree(buf_gpu);
    cudaFree(norm_gpu);

    return k;
}
