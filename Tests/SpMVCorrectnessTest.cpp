#include "Solver/Kernels/mathfunc.h"

int main() {
    std::vector<int> ia = {0, 5, 10, 15, 20, 25};
    std::vector<int> ja = {
        0, 1, 2, 3, 4,
        0, 1, 2, 3, 4,
        0, 1, 2, 3, 4,
        0, 1, 2, 3, 4,
        0, 1, 2, 3, 4};
    std::vector<double> a = {
        11, 1, 2, 3, 4,
        1, 12, 1, 2, 3,
        2, 1, 13, 1, 2,
        3, 2, 1, 14, 1,
        4, 3, 2, 1, 15};

    std::vector<double> x = {1, 2, 3, 4, 5};
    std::vector<double> res(x.size());
    double norm_res, norm_ans;

    // SpMV solution
    spMV(ia.data(), ja.data(), a.data(), x.data(), ia.size(), res.data());
    dot(res.data(), res.data(), res.size(), norm_res);

    //Correct answer
    std::vector<double> ans = {51, 51, 57, 71, 95};
    dot(ans.data(), ans.data(), ans.size(), norm_ans);

    return norm_res != norm_ans;
}
