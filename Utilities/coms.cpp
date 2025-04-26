#include "coms.h"

#include <algorithm>
#include <cassert>
#include <mpi.h>
#include <unordered_map>
#include <vector>

// #include "input.h"


void ComInit(const std::vector<int> &ia, const std::vector<int> &ja, const std::vector<int> &Part, const std::vector<int> &L2G, const int
             MyID, std::vector<int> &RecvOffset, std::vector<int> &SendOffset, std::vector<int> &Send, std::vector<int> &Recv, std::vector<int> &Neighbors) {

    int NumProc;
    MPI_Comm_size(MPI_COMM_WORLD,&NumProc);

    std::vector<std::vector<int>> SendToProcess(NumProc);
    std::vector<std::vector<int>> RecvFromProcess(NumProc);

    std::unordered_map<int, int> P2N;

    for (unsigned int i = 0; i < ia.size() - 1; i++) {
        for (int col = ia[i]; col < ia[i + 1]; col++) {
            const int j = ja[col];
            if (Part[j] != MyID) {
                if (P2N.find(Part[j]) == P2N.end()) {
                    P2N[Part[j]] = P2N.size();
                }

                SendToProcess[P2N[Part[j]]].push_back(i);
                RecvFromProcess[P2N[Part[j]]].push_back(j);
            }
        }
    }

    SendOffset.push_back(0);
    RecvOffset.push_back(0);

    Neighbors.resize(P2N.size());
    for (const auto p : P2N) {
        Neighbors[p.second] = p.first;
    }

    for (unsigned int p = 0; p < Neighbors.size(); p++) {
        std::sort(SendToProcess[p].begin(), SendToProcess[p].end(),
            [&L2G](const int i, const int j) {
                return L2G[i] < L2G[j];
            });

        int k = 1;
        Send.push_back(SendToProcess[p][0]);
        for (unsigned int i = 1; i < SendToProcess[p].size(); i++) {
            if (SendToProcess[p][i] != SendToProcess[p][i - 1]) {
                Send.push_back(SendToProcess[p][i]);
                k++;
            }
        }
        SendOffset.push_back(SendOffset.back() + k);
    }

    for (unsigned int p = 0; p < Neighbors.size(); p++) {
        std::sort(RecvFromProcess[p].begin(), RecvFromProcess[p].end(),
        [&L2G](const int i, const int j) {
            return L2G[i] < L2G[j];
        });

        int k = 1;
        Recv.push_back(RecvFromProcess[p][0]);
        for (unsigned int i = 1; i < RecvFromProcess[p].size(); i++) {
            if (RecvFromProcess[p][i] != RecvFromProcess[p][i - 1]) {
                Recv.push_back(RecvFromProcess[p][i]);
                k++;
            }
        }
        RecvOffset.push_back(RecvOffset.back() + k);
    }
}


void ComUpdate(double *b, const std::vector<int> &RecvOffset, const std::vector<int> &SendOffset,
               const std::vector<int> &Neighbors, const std::vector<int> &Send, const
               std::vector<int> &Recv) {

    std::vector<double> SendBuf, RecvBuf;
    SendBuf.resize(Send.size());
    RecvBuf.resize(Recv.size());

    assert(static_cast<int>(Recv.size()) == RecvOffset.back() && "Recv size incorrect");
    assert(static_cast<int>(Send.size()) == SendOffset.back() && "Send size incorrect");

    std::vector<MPI_Request> Request;
    std::vector<MPI_Status> Status;
    Request.resize(2 * Neighbors.size());
    Status.resize(2 * Neighbors.size());

    int nreq = 0;

    for (unsigned int p = 0; p < Neighbors.size(); p++) {
        const int size = RecvOffset[p + 1] - RecvOffset[p];
        const int neighbor_id = Neighbors[p];
        const int mpires = MPI_Irecv(&RecvBuf[RecvOffset[p]], size, MPI_DOUBLE, neighbor_id, 0, MPI_COMM_WORLD, &(Request[nreq]));
        assert(mpires == MPI_SUCCESS && "MPI_Irecv failed");
        nreq++;
    }

    for (unsigned int i = 0; i < Send.size(); i++) {
        SendBuf[i] = b[Send[i]];
    }

    for (unsigned int p = 0; p < Neighbors.size(); p++) {
        const int size = SendOffset[p + 1] - SendOffset[p];
        const int neighbor_id = Neighbors[p];
        const int mpires = MPI_Isend(&SendBuf[SendOffset[p]], size, MPI_DOUBLE, neighbor_id, 0, MPI_COMM_WORLD, &(Request[nreq]));
        assert(mpires == MPI_SUCCESS && "MPI_Isend failed");
        nreq++;
    }

    if (nreq > 0) {
        const int mpires = MPI_Waitall(nreq, &Request[0], &Status[0]);
        assert(mpires == MPI_SUCCESS && "MPI_Waitall failed");
    }

    for (unsigned int i = 0; i < Recv.size(); i++) {
        b[Recv[i]] = RecvBuf[i];
    }
}
