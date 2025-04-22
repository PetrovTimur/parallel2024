if (USE_MPI)
    add_library(Coms SHARED coms.cpp coms.h)
endif ()

add_library(IOs SHARED
        input.cpp
        input.h
        argparse.h
        logger.h
        logger.cpp)