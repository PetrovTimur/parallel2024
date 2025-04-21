#include <iostream>

int main() {
    cudaDeviceProp prop{};
    int device;
    cudaGetDevice(&device);
    cudaGetDeviceProperties(&prop, device);

    std::cout << "Device " << prop.name << std::endl;
    std::cout << prop.major << "." << prop.minor << std::endl;
    std::cout << prop.maxThreadsPerBlock << std::endl;
    std::cout << prop.maxBlocksPerMultiProcessor << std::endl;
    std::cout << prop.multiProcessorCount <<std::endl;
    std::cout << prop.maxThreadsPerMultiProcessor << std::endl;

    std::cout << std::endl << "-----------------" << std::endl << std::endl;

}
