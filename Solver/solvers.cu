#include "cuda_runtime.h"
#include "thrust/device_vector.h"
#include "Solver/Kernels/mathfunc.cuh"
#include <Utilities/logger.h>
#include <Utilities/cuda_helper.cuh>
#include "solvers.cuh"
#include "omp.h"


int solve(const int *ia, const int *ja, const float *a, const float *b, const float *diag, int size, float *res, const double eps, const int maxit) {
    const int N = size;
    int blocks = 204;
    int threads = 256;

    getDeviceSpecs(blocks, threads);

    LOG_INFO << "Using " << blocks << " blocks, " << threads << " threads/block" << std::endl;

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

    // CUDA timing events
    cudaEvent_t start_solve, end_solve, start_iter, end_iter;
    cudaEventCreate(&start_solve);
    cudaEventCreate(&end_solve);
    cudaEventCreate(&start_iter);
    cudaEventCreate(&end_iter);

    cudaEventRecord(start_solve);

    do {
        k++;

        LOG_INFO << "Iteration " << k << std::endl;

        cudaEventRecord(start_iter);

        multiply<<<blocks, threads>>>(r.data().get(), diag_gpu.data().get(), z.data().get(), N);
        getLastCudaError("Kernel execution failed");

        rho[0] = rho[1];

        assert(r.size() == z.size() && r.size() == N && "Size mismatch");
        dot_gpu(threads, blocks, r.data().get(), z.data().get(), vec_buf.data().get(), rho.data().get() + 1, N);
        getLastCudaError("Kernel execution failed");

        if (k == 1)
            p = z;
        else {
            const float beta = rho[1] / rho[0];
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

        cudaEventRecord(end_iter);
        cudaEventSynchronize(end_iter);
        float iter_ms = 0.0f;
        cudaEventElapsedTime(&iter_ms, start_iter, end_iter);
        LOG_INFO << "Time " << iter_ms / 1000.0 << std::endl;

        dot_gpu(threads, blocks, x.data().get(), x.data().get(), vec_buf.data().get(), norm_gpu, N);
        cudaMemcpy(&norm, norm_gpu, sizeof(float), cudaMemcpyDeviceToHost);
        LOG_INFO << "Solution norm: " << std::sqrt(norm) << std::endl;
        LOG << "--------------------------------------------------" << std::endl << std::endl;

    }
    while (rho[1] > eps * eps && k < maxit);

    cudaEventRecord(end_solve);
    cudaEventSynchronize(end_solve);
    float solve_ms = 0.0f;
    cudaEventElapsedTime(&solve_ms, start_solve, end_solve);
    LOG_INFO << "Solve time " << solve_ms / 1000.0 << std::endl;

    cudaMemcpy(res, x.data().get(), N * sizeof(float), cudaMemcpyDeviceToHost);

    cudaFree(buf_gpu);
    cudaFree(norm_gpu);

    cudaEventDestroy(start_solve);
    cudaEventDestroy(end_solve);
    cudaEventDestroy(start_iter);
    cudaEventDestroy(end_iter);

    return k;
}
