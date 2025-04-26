#include "omp.h"
#include "csr.h"
#include "Utilities/input.h"
#include "Utilities/logger.h"
#include "Utilities/argparse.h"
#include "solvers.h"
#include <mpi.h>
#include <unordered_map>


int main(int argc, char** argv) {
    int mpi_res = MPI_Init(&argc, &argv); // первым делом подключаемся к MPI
    if(mpi_res != MPI_SUCCESS) {
        MPI_Abort(MPI_COMM_WORLD, -1);
        exit(mpi_res);
    }

    int NumProc;
    mpi_res = MPI_Comm_size(MPI_COMM_WORLD,&NumProc); // узнаем число процессов
    if(mpi_res != MPI_SUCCESS){
        MPI_Abort(MPI_COMM_WORLD, -1);
        exit(mpi_res);
    }

    int MyID;
    mpi_res = MPI_Comm_rank(MPI_COMM_WORLD,&MyID); // узнаем номер данного процесса
    if(mpi_res != MPI_SUCCESS){
        MPI_Abort(MPI_COMM_WORLD, -1);
        exit(mpi_res);
    }

    arguments args{};

    args.log_dir = "";
    args.eps = 1e-3;
    args.maxit = 100;

    argp_parse (&argp, argc, argv, 0, nullptr, &args);

    if (!args.log_dir.empty() && (MyID == 0)) {
        Logger::setLogDirectory(args.log_dir);
        LOG_INFO << "Starting with log directory: " << args.log_dir << std::endl;
    }

    omp_set_num_threads(1);
    // omp_set_num_threads(omp_get_max_threads());

    int Nx = args.Nx;
    checkInput(Nx, "Nx");
    int Ny = args.Ny;
    checkInput(Ny, "Ny");
    int K1 = args.K1;
    checkInput(K1, "K1");
    int K2 = args.K2;
    checkInput(K2, "K2");
    int Px = args.Px;
    checkInput(Px, "Px");
    int Py = args.Py;
    checkInput(Py, "Py");

    double eps = args.eps;
    checkInput(eps, "eps");
    int maxit = args.maxit;
    checkInput(maxit, "maxit");


    if (MyID == 0) {
        LOG_INFO << "Nx = " << Nx << ", Ny = " << Ny << ", K1 = " << K1 << ", K2 = " << K2 << ", Px = " << Px << ", Py = " << Py << std::endl;
    }

    if (NumProc != Px * Py) {
        if (MyID == 0) {
            LOG_ERROR << "Number of processes Px*Py doesn't match MPI config" << std::endl;
        }

        MPI_Finalize();
        return 0;
    }

    // std::cout << omp_get_max_threads() << std::endl;
    // omp_set_num_threads(omp_get_max_threads() / P0);
    // std::cout << omp_get_max_threads() << std::endl;

    std::vector<int> L2G;
    // std::vector<int> G2L;
    std::vector<int> Part;

    input(Nx, Ny, Px, Py, MyID, L2G, Part);


    // std::vector<int> ia(N0 + 1);
    std::vector<int> ia_en;
    // std::vector<int> ja(7 * (N0 + 1));
    std::vector<int> ja_en;

    // std::cout << "My ID = " << MyID << std::endl;
    // std::cout << "L2G size = " << L2G.size() << std::endl;
    // printVector(L2G, std::cout);
    //
    // MPI_Barrier(MPI_COMM_WORLD);

    double start = MPI_Wtime();
    makeIncidenceMatrixCSR(Nx, Ny, K1, K2, L2G, ia_en, ja_en, Part);
    MPI_Barrier(MPI_COMM_WORLD);
    double end = MPI_Wtime();

    if (MyID == 0) {
        LOG_INFO << "EN construction done in " << end - start << " seconds" << std::endl;
    }


    std::unordered_map<int, int> G2L;
    fillG2L(L2G, G2L);

    std::unordered_map<int, int> G2L_nodes;
    constructG2L(ia_en, ja_en, G2L_nodes);

    localizeCSR(ia_en.data(), ia_en.size(), ja_en.data(), G2L_nodes);

    // printVector(L2G, std::cout);
    // std::cout << "G2L_nodes size = " << G2L_nodes.size() << std::endl;
    // printVector(ja_en, std::cout);

    // std::vector<int> ia_ne, ja_ne;
    int *ia_ne, *ja_ne;

    start = MPI_Wtime();
    transposeCSR(ia_en, ja_en, G2L_nodes.size(), ia_ne, ja_ne);
    MPI_Barrier(MPI_COMM_WORLD);
    end = MPI_Wtime();

    if (MyID == 0) {
        LOG_INFO << "Transpose (NE) construction done in " << end - start << " seconds" << std::endl;
    }

    // printVector(ja_ne, std::cout);

    // MPI_Barrier(MPI_COMM_WORLD);
    // std::cout << "My ID = " << MyID << ", Part size " << Part.size() << ", Ne = " << ia_en.size() - 1 << std::endl;
    // printVector(Part, std::cout);

    std::vector<int> ia_ee, ja_ee;

    start = MPI_Wtime();
    buildAdjacencyMatrixCSRUsingSort(ia_en.data(), ja_en.data(), ia_ne, ja_ne, ia_ee, ja_ee, ia_en.size() - 1, Part, MyID);
    MPI_Barrier(MPI_COMM_WORLD);
    end = MPI_Wtime();

    if (MyID == 0) {
        LOG_INFO << "Adjacency matrix (EE) construction done in " << end - start << " seconds" << std::endl;
    }

    std::vector<double> a(ja_ee.size());
    std::vector<double> b(ia_ee.size() - 1);
    std::vector<double> diag(ia_ee.size() - 1);

    start = MPI_Wtime();
    fillCSR(ia_ee, ja_ee, L2G, a, b, diag);
    MPI_Barrier(MPI_COMM_WORLD);
    end = MPI_Wtime();

    if (MyID == 0) {
        LOG_INFO << "Fill done in " << end - start << " seconds" << std::endl << std::endl;
    }


    // std::cout << "My ID = " << MyID << std::endl;
    // printVector(ja_ee, std::cout);

    // printVector(ia);
    // printVector(ja);
    // if (MyID == 0) {
    //     LOG_INFO << "My ID = " << MyID << ", b size: " << b.size() << std::endl;
    //     printVector(ia_ee, LOG);
    //     printVector(a, LOG);
    //     printVector(b, LOG);
    // }
    // printVector(b);

    std::vector<double> res(ia_ee.size() - 1);

    start = MPI_Wtime();
    int iterations = solve(MyID, Part, L2G, ia_ee, ja_ee, a, b, diag, res);
    MPI_Barrier(MPI_COMM_WORLD);
    end = MPI_Wtime();

    if (MyID == 0) {
        LOG_INFO << "Work took " << end - start << " seconds" << std::endl;
        LOG_INFO << "Convergence required "  << iterations << " iterations" << std::endl;
        LOG_INFO << res[0] << std::endl;
    }

    delete[] ia_ne;
    delete[] ja_ne;


    MPI_Finalize();

    return 0;
}
