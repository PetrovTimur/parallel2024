cmake_minimum_required(VERSION 3.12)
project(parallel2024_task1)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin" CACHE INTERNAL "Where to place executables")

option(ENABLE_TESTS "Enable automatic testing" ON)
if (ENABLE_TESTS)
    enable_testing()
    add_subdirectory(Tests)
endif ()


option(USE_DEBUG_MODE "Use debug tools" OFF)
if (USE_DEBUG_MODE)
    message("LULW")
    add_definitions(-DUSE_DEBUG_MODE)
endif ()


set(TOPLEVEL_INCLUDE_DIR ${PROJECT_SOURCE_DIR})

add_subdirectory(Utilities)
add_subdirectory(Solver)
