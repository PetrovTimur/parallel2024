add_subdirectory(Kernels)

if (USE_CUDA)
    add_executable(CGSolver CGSolver.cpp
            csr.cpp
            csr.h
            solvers.cu
            solvers.cuh)

    add_dependencies(CGSolver IOs MathFuncCUDA)
    target_link_libraries(CGSolver PUBLIC IOs MathFuncCUDA)
else ()
    if (USE_MPI)
        set(source CGSolverMPI.cpp)
    else ()
        set(source CGSolver.cpp)
    endif ()

    add_executable(CGSolver ${source}
            csr.cpp
            csr.h
            solvers.cpp
            solvers.h)

    add_dependencies(CGSolver IOs MathFunc)
    target_link_libraries(CGSolver PUBLIC IOs MathFunc)
    if (USE_MPI)
        add_dependencies(CGSolver Coms)
        target_link_libraries(CGSolver PUBLIC Coms)
    endif ()
endif()