#include "Solver/Kernels/mathfunc.h"
#include <iostream>
#include <vector>

bool Test1() {
    std::vector<int> arr = {0, 2, 3, 5, 5, 2, 6, 7, 2, 7, 2, 7, 9, 4, 12, 0, 0, 345, 234, 12};  // Row pointers
    std::vector<int> expected_ans(arr.size(), arr[0]);

    for (unsigned int i = 1; i < arr.size(); i++) {
        expected_ans[i] = expected_ans[i - 1] + arr[i];
    }

    scan(arr.data(), arr.data(), arr.size());

    bool success = true;

    // Check ans values
    for (unsigned int i = 0; i < arr.size(); i++) {
        if (arr[i] != expected_ans[i]) {
            std::cout << "Error in ans[" << i << "]: expected " << expected_ans[i]
                      << ", got " << arr[i] << std::endl;
            success = false;
        }
    }

    return success;
}


int main() {
    bool all_passed = true;

    std::cout << "Testing parallel scan... Test 1 ";
    all_passed &= Test1();
    std::cout << (all_passed ? "PASSED\n" : "FAILED\n");

    return all_passed ? 0 : 1;
}