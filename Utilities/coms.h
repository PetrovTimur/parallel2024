#ifndef COMS_H
#define COMS_H
#include <mpi.h>
#include <vector>

void ComInitOffsets(int top_halo, int left_halo, int right_halo, int bottom_halo, int i_count, int j_count, std::vector<int> &recv_offset, std::vector<int> &send_offset);

void Com(int MyID, int Px, int top_halo, int left_halo, int right_halo, int bottom_halo, int i_count, int j_count,
         const double *b, std::vector<int> &recv_offset, std::vector<int> &send_offset,
         std::vector<double> &recv_buf, std::vector<double> &send_buf);

#endif //COMS_H
