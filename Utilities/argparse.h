#ifndef ARGPARSE_H
#define ARGPARSE_H

#include <stdlib.h>
#include <argp.h>

const char *argp_program_version =
  "CGSolver 1.0";
const char *argp_program_bug_address =
  "<s02240520@stud.cs.msu.ru>";

/* Program documentation. */
static char doc[] =
  "CGSolver 20 20 5 6 -- a program with options and arguments";

/* A description of the arguments we accept. */
static char args_doc[] = "Nx Ny K1 K2";

/* The options we understand. */
static struct argp_option options[] = {
    // {"verbose",  'v', 0,      0,  "Produce verbose output" },
    // {"quiet",    'q', 0,      0,  "Don't produce any output" },
    // {"silent",   's', 0,      OPTION_ALIAS },
    {"output",   'o', "FILE", 0,
     "Output to FILE instead of standard output" },
    { 0 }
};

/* Used by main to communicate with parse_opt. */
struct arguments
{
    char *args[4];                /* arg1 & arg2 */
    const char *output_file;
};

/* Parse a single option. */
static error_t parse_opt(int key, char *arg, struct argp_state *state) {
    /* Get the input argument from argp_parse, which we
       know is a pointer to our arguments structure. */
    struct arguments *arguments = (struct arguments *) state->input;

    switch (key)
    {
        // case 'q': case 's':
        //     arguments->silent = 1;
        // break;
        case 'o':
            arguments->output_file = arg;
        break;

        case ARGP_KEY_ARG:
            if (state->arg_num > 4)
                /* Too many arguments. */
                    argp_usage (state);

        arguments->args[state->arg_num] = arg;

        break;

        case ARGP_KEY_END:
            if (state->arg_num < 4)
                /* Not enough arguments. */
                    argp_usage (state);
        break;

        default:
            return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

/* Our argp parser. */
static struct argp argp = { options, parse_opt, args_doc, doc };

#endif //ARGPARSE_H
