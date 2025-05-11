#ifndef CUDAUTILS_CUH
#define CUDAUTILS_CUH

#include <iostream>

#define checkCudaErrors(val) check((val), #val, __FILE__, __LINE__)

inline void check(cudaError_t err, const char* const func, const char* const file,
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

inline void getDeviceSpecs(int &blocks, int &threads) {
    cudaDeviceProp prop{};
    int device;
    cudaGetDevice(&device);
    cudaGetDeviceProperties(&prop, device);

    #ifdef DEBUG_MODE
    LOG_DEBUG << "Device " << prop.name << ", " << prop.major << "." << prop.minor << std::endl;
    LOG_DEBUG << "Max threads per block: " << prop.maxThreadsPerBlock << std::endl;
    LOG_DEBUG << "Max blocks per SM: " << prop.maxBlocksPerMultiProcessor << std::endl;
    LOG_DEBUG << "SM count: " << prop.multiProcessorCount << std::endl;
    LOG_DEBUG << "Max threads per SM: " << prop.maxThreadsPerMultiProcessor << std::endl;

    LOG_DEBUG << "-----------------" << std::endl << std::endl;
    #endif

    blocks = prop.multiProcessorCount * prop.maxThreadsPerMultiProcessor / 256;
    threads = 256;
}


#endif //CUDAUTILS_CUH
