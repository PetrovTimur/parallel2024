cmake_minimum_required(VERSION 3.29)
project(parallel2024_task1)

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

add_executable(parallel2024_task1 main.cpp)
