#ifndef SOLVERS_H
#define SOLVERS_H
#ifdef USE_MPI
#include <mpi.h>
#endif
#include <vector>

#ifdef USE_MPI
int solve(int MyID, int Px, int top_halo, int left_halo, int right_halo, int bottom_halo, int i_count, int j_count,
        std::vector<int> &recv_offset, std::vector<int> &send_offset,
        std::vector<double> &recv_buf, std::vector<double> &send_buf, std::vector<MPI_Request> &recv_req,
        std::vector<MPI_Request> &send_req, std::vector<MPI_Status> &recv_stat, std::vector<MPI_Status> &send_stat,
        std::vector<int> &ia, std::vector<int> &ja, std::vector<double> &a, std::vector<double> &b,
        std::vector<double> &diag, std::vector<double> &res);
#else
int solve(std::vector<int> &ia, std::vector<int> &ja, std::vector<double> &a, std::vector<double> &b,
          std::vector<double> &diag, std::vector<double> &res);
#endif


#endif //SOLVERS_H
