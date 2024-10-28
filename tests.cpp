#include <iostream>
#include <ostream>

#include "mathfunc.cpp"


int main() {
    std::vector<double> x = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
    std::vector<double> y = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};

    double res = dot(x, y);

    std::cout << res << std::endl;

    return 0;
}