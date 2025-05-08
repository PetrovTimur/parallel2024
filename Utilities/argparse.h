#ifndef ARGPARSE_H
#define ARGPARSE_H

#include <cstdlib>
#include <argp.h>
#include <string>

const char * argp_program_version = "CGSolver 1.0";
const char * argp_program_bug_address = "<s02240520@stud.cs.msu.ru>";

/* Program documentation. */
#if defined(USE_MPI)
static char doc[] = "CGSolver -- a parallel solver using MPI";
#elif defined(USE_CUDA)
static char doc[] = "CGSolver -- a parallel solver using CUDA";
#else
static char doc[] = "CGSolver -- a parallel solver using OpenMP";
#endif


/* A description of the arguments we accept. */
#ifdef USE_MPI
static char args_doc[] = "[--Nx VALUE] [--Ny VALUE] [--K1 VALUE] [--K2 VALUE] [--Px VALUE] [--Py VALUE]";
#else
static char args_doc[] = "[--Nx VALUE] [--Ny VALUE] [--K1 VALUE] [--K2 VALUE]";
#endif

/* The options we understand. */
static struct argp_option options[] = {
    {"log_dir", 'd', "DIR", 0, "Output to DIR instead of standard log directory"},
    {"eps", 'e', "VALUE", 0, "Precision of calculations"},
    {"maxit", 'i', "VALUE", 0, "Maximum iterations"},
    {"Nx", 1001, "VALUE", 0, "X dimension"},
    {"Ny", 1002, "VALUE", 0, "Y dimension"},
    {"K1", 1003, "VALUE", 0, "K1 parameter"},
    {"K2", 1004, "VALUE", 0, "K2 parameter"},
#ifdef USE_MPI
    {"Px", 1005, "VALUE", 0, "X processes"},
    {"Py", 1006, "VALUE", 0, "Y processes"},
#endif
    {0}
};

/* Used by main to communicate with parse_opt. */
struct arguments {
    int Nx, Ny, K1, K2;
#ifdef USE_MPI
    int Px, Py;
#endif
    float eps;
    int maxit;
    std::string log_dir;
};

/* Parse a single option. */
static error_t parse_opt(int key, char *arg, struct argp_state *state) {
    /* Get the input argument from argp_parse, which we
       know is a pointer to our arguments structure. */
    auto *arguments = static_cast<struct arguments *>(state->input);

    switch (key)
    {
        // case 'q': case 's':
        //     arguments->silent = 1;
        // break;
        case 'd':
            arguments->log_dir = arg;
        break;

        case 'e':
            arguments->eps = std::stof(arg);
        break;

        case 'i':
            arguments->maxit = std::stoi(arg);
        break;

        case 1001:  // Nx
            arguments->Nx = std::stoi(arg);
        break;
        case 1002:  // Ny
            arguments->Ny = std::stoi(arg);
        break;
        case 1003:  // K1
            arguments->K1 = std::stoi(arg);
        break;
        case 1004:  // K2
            arguments->K2 = std::stoi(arg);
        break;
        #ifdef USE_MPI
        case 1005:  // Px
            arguments->Px = std::stoi(arg);
        break;
        case 1006:  // Py
            arguments->Py = std::stoi(arg);
        break;
        #endif

        case ARGP_KEY_END:
            // Check that all required arguments were provided
            if (arguments->Nx == 0 || arguments->Ny == 0 ||
                arguments->K1 == 0 || arguments->K2 == 0
                #ifdef USE_MPI
                || arguments->Px == 0 || arguments->Py == 0
                #endif
            ) {
                argp_usage(state);
            }
        break;
        default:
            return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

/* Our argp parser. */
static struct argp argp = { options, parse_opt, args_doc, doc };

#endif //ARGPARSE_H
