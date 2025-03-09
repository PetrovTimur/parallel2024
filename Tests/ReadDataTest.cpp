#include "Utilities/input.h"
#include <iostream>
#include <cassert>

void testReadData() {
    try {
        auto result = readData(PROJECT_SOURCE_DIR "/Data/mesh100K/elements.txt", PROJECT_SOURCE_DIR "/Data/mesh100K/elements.dat");
        int Nn = std::get<0>(result);
        int Ne = std::get<1>(result);
        int M = std::get<2>(result);
        int* data = std::get<3>(result);

        // Check that the data is not null
        assert(data != nullptr);

        // Clean up allocated memory
        delete[] data;

        std::cout << "Test passed: No errors thrown and data is valid." << std::endl;
    } catch (const std::exception &e) {
        std::cerr << "Test failed: Exception thrown - " << e.what() << std::endl;
    }
}

int main() {
    testReadData();
    return 0;
}