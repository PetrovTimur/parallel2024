add_subdirectory(Kernels)

add_executable(Solver main.cpp
        csr.cpp
        csr.h
        solvers.cpp
        solvers.h)

add_dependencies(Solver Utilities MathFunc)
target_link_libraries(Solver PUBLIC MathFunc Utilities)