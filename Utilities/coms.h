#ifndef COMS_H
#define COMS_H
#include <mpi.h>
#include <vector>

void ComInit(const std::vector<int> &ia, const std::vector<int> &ja, const std::vector<int> &Part, const std::vector<int> &L2G, int
                    MyID, std::vector<int> &RecvOffset, std::vector<int> &SendOffset, std::vector<int> &Send, std::vector<int> &Recv, std
                    ::vector<int> &Neighbors);

void ComUpdate(
    double *b, const std::vector<int> &RecvOffset, const std::vector<int> &SendOffset,
    const std::vector<int> &Neighbors, const std::vector<int> &Send, const
    std::vector<int> &Recv);

#endif //COMS_H
