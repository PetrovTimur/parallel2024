add_subdirectory(Kernels)

add_executable(CGSolver main.cpp
        csr.cpp
        csr.h
        solvers.cpp
        solvers.h)

add_dependencies(CGSolver Utilities MathFunc)
target_link_libraries(CGSolver PUBLIC MathFunc Utilities)