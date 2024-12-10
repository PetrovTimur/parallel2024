#include "coms.h"
#include <mpi.h>
#include <vector>


void ComInitOffsets(int top_halo, int left_halo, int right_halo, int bottom_halo, int i_count, int j_count, std::vector<int> &recv_offset, std::vector<int> &send_offset) {
    recv_offset[0] = 0;
    recv_offset[1] = top_halo * j_count;
    recv_offset[2] = recv_offset[1] + top_halo * right_halo;
    recv_offset[3] = recv_offset[2] + left_halo * i_count;
    recv_offset[4] = recv_offset[3] + right_halo * i_count;
    recv_offset[5] = recv_offset[4] + left_halo * bottom_halo;
    recv_offset[6] = recv_offset[5] + bottom_halo * j_count;

    send_offset[0] = 0;
    send_offset[1] = top_halo * j_count;
    send_offset[2] = send_offset[1] + top_halo * right_halo;
    send_offset[3] = send_offset[2] + left_halo * i_count;
    send_offset[4] = send_offset[3] + right_halo * i_count;
    send_offset[5] = send_offset[4] + left_halo * bottom_halo;
    send_offset[6] = send_offset[5] + bottom_halo * j_count;
}


void Com(int MyID, int Px, int top_halo, int left_halo, int right_halo, int bottom_halo, int i_count, int j_count,
        std::vector<double> &b, std::vector<int> &recv_offset, std::vector<int> &send_offset,
        std::vector<double> &recv_buf, std::vector<double> &send_buf) {
    int nrreq = 0;
    int mpires;

    std::vector<MPI_Request> recv_req(top_halo + top_halo*right_halo + left_halo + right_halo + left_halo*bottom_halo + bottom_halo);
    std::vector<MPI_Status> recv_stat(top_halo + top_halo*right_halo + left_halo + right_halo + left_halo*bottom_halo + bottom_halo);
    std::vector<MPI_Request> send_req(top_halo + top_halo*right_halo + left_halo + right_halo + left_halo*bottom_halo + bottom_halo);
    std::vector<MPI_Status> send_stat(top_halo + top_halo*right_halo + left_halo + right_halo + left_halo*bottom_halo + bottom_halo);

    if (top_halo) {
        int size = recv_offset[1]-recv_offset[0];
        MPI_Irecv(&recv_buf[recv_offset[0]], size, MPI_DOUBLE, MyID - Px, 0, MPI_COMM_WORLD, &recv_req[nrreq]);
        nrreq++;

        if (right_halo) {
            int size = recv_offset[2]-recv_offset[1];
            MPI_Irecv(&recv_buf[recv_offset[1]], size, MPI_DOUBLE, MyID - Px + 1, 0, MPI_COMM_WORLD, &recv_req[nrreq]);
            nrreq++;
        }
    }

    if (left_halo) {
        int size = recv_offset[3]-recv_offset[2];
        MPI_Irecv(&recv_buf[recv_offset[2]], size, MPI_DOUBLE, MyID - 1, 0, MPI_COMM_WORLD, &recv_req[nrreq]);
        nrreq++;
    }

    if (right_halo) {
        int size = recv_offset[4]-recv_offset[3];
        MPI_Irecv(&recv_buf[recv_offset[3]], size, MPI_DOUBLE, MyID + 1, 0, MPI_COMM_WORLD, &recv_req[nrreq]);
        nrreq++;
    }

    if (bottom_halo) {
        if (left_halo) {
            int size = recv_offset[5]-recv_offset[4];
            MPI_Irecv(&recv_buf[recv_offset[4]], size, MPI_DOUBLE, MyID + Px - 1, 0, MPI_COMM_WORLD, &recv_req[nrreq]);
            nrreq++;
        }

        int size = recv_offset[6]-recv_offset[5];
        MPI_Irecv(&recv_buf[recv_offset[5]], size, MPI_DOUBLE, MyID + Px, 0, MPI_COMM_WORLD, &recv_req[nrreq]);
        nrreq++;
    }

    int p = 0;

    int nsreq = 0;
    if (top_halo) {
        #pragma omp parallel for proc_bind(spread)
        for (int j = 0; j < j_count; j++) {
            send_buf[p] = b[j];
            p++;
        }

        int size = send_offset[1]-send_offset[0];
        MPI_Isend(&send_buf[send_offset[0]], size, MPI_DOUBLE, MyID - Px, 0, MPI_COMM_WORLD, &send_req[nsreq]);
        nsreq++;

        if (right_halo) {
            send_buf[p] = b[j_count - 1];
            p++;

            int size = send_offset[2]-send_offset[1];
            MPI_Isend(&send_buf[send_offset[1]], size, MPI_DOUBLE, MyID - Px + 1, 0, MPI_COMM_WORLD, &send_req[nsreq]);
            nsreq++;
        }
    }

    if (left_halo) {
        #pragma omp parallel for proc_bind(spread)
        for (int i = 0; i < i_count; i++) {
            send_buf[p] = b[i * j_count];
            p++;
        }

        int size = send_offset[3]-send_offset[2];
        MPI_Isend(&send_buf[send_offset[2]], size, MPI_DOUBLE, MyID - 1, 0, MPI_COMM_WORLD, &send_req[nsreq]);
        nsreq++;
    }

    if (right_halo) {
        #pragma omp parallel for proc_bind(spread)
        for (int i = 0; i < i_count; i++) {
            send_buf[p] = b[i * j_count + j_count - 1];
            p++;
        }

        int size = send_offset[4]-send_offset[3];
        MPI_Isend(&send_buf[send_offset[3]], size, MPI_DOUBLE, MyID + 1, 0, MPI_COMM_WORLD, &send_req[nsreq]);
        nsreq++;
    }

    if (bottom_halo) {
        if (left_halo) {
            send_buf[p] = b[(i_count - 1) * j_count];
            p++;

            int size = send_offset[5]-send_offset[4];
            MPI_Isend(&send_buf[send_offset[4]], size, MPI_DOUBLE, MyID + Px - 1, 0, MPI_COMM_WORLD, &send_req[nsreq]);
            nsreq++;
        }

        #pragma omp parallel for proc_bind(spread)
        for (int j = 0; j < j_count; j++) {
            send_buf[p] = b[(i_count - 1) * j_count + j];
            p++;
        }

        int size = send_offset[6]-send_offset[5];
        MPI_Isend(&send_buf[send_offset[5]], size, MPI_DOUBLE, MyID + Px, 0, MPI_COMM_WORLD, &send_req[nsreq]);
        nsreq++;
    }

    // if (MyID == 9) {
    //     std::cout << "Offsets for send" << std::endl;
    //     for (int i : send_offset)
    //         std::cout << i << " ";
    //     std::cout << std::endl;
    // }



    if (nrreq>0){ // ждем завершения получения
        mpires = MPI_Waitall(nrreq, &recv_req[0], &recv_stat[0]);
        // ASSERT(mpires==MPI_SUCCESS, "MPI_Waitall (recv) failed");
        // std::cout << (mpires == MPI_SUCCESS) << std::endl;
    }

    if (nsreq>0){ // ждем завершения отправок
        mpires = MPI_Waitall(nsreq, &send_req[0], &send_stat[0]);
        // ASSERT(mpires==MPI_SUCCESS, "MPI_Waitall (recv) failed");
        // std::cout << (mpires == MPI_SUCCESS) << std::endl;
    }

    // std::cout << "After recv: MyID: " << MyID << std::endl;
    // std::cout << "rec buf size: " << recv_buf.size() << std::endl;
    // for (double x : recv_buf) {
    //     std::cout << x << std::endl;
    // }
    // std::cout << std::endl;
}
