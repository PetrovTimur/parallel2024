cmake_minimum_required(VERSION 3.29)
project(parallel2024_task1)

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")


add_executable(parallel2024_task1 main.cpp
        mathfunc.cpp
        mathfunc.h
        csr.cpp
        csr.h
        input.cpp
        input.h
        solvers.cpp
        solvers.h)


#testing binary
add_executable(math_test tests.cpp)
#target_link_libraries(cpp_test PRIVATE sum_integers)

enable_testing()

add_test(
        NAME math_test
        COMMAND $<TARGET_FILE:math_test>
)