cmake_minimum_required(VERSION 3.29)
project(parallel2024_task1)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
enable_testing()

set(TOPLEVEL_INCLUDE_DIR ${PROJECT_SOURCE_DIR})

add_executable(parallel2024_task1 main.cpp
        csr.cpp
        csr.h
        input.cpp
        input.h
        solvers.cpp
        solvers.h)

add_subdirectory(Tests)
add_subdirectory(Ops)

target_link_libraries(parallel2024_task1 MathFunc)
