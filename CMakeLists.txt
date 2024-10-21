cmake_minimum_required(VERSION 3.29)
project(parallel2024_task1)

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")


#add_library(mathfunc mathfunc.cpp)

add_executable(parallel2024_task1 main.cpp
        mathfunc.cpp
        mathfunc.h
        csr.cpp
        csr.h
        input.cpp
        input.h
        solvers.cpp
        solvers.h)

#target_link_libraries(parallel2024_task1 mathfunc)
#add_dependencies(parallel2024_task1 )