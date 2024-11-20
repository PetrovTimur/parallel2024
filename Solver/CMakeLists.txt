add_subdirectory(Kernels)

add_executable(CGSolver main.cpp
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