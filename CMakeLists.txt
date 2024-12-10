cmake_minimum_required(VERSION 3.12)
project(parallel2024_task1)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin" CACHE INTERNAL "Where to place executables")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/lib" CACHE INTERNAL "Where to place libraries")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/lib" CACHE INTERNAL "Where to place libraries")

option(ENABLE_TESTS "Enable automatic testing" ON)


option(USE_DEBUG_MODE "Use debug tools" OFF)
if (USE_DEBUG_MODE)
    add_definitions(-DUSE_DEBUG_MODE)
endif ()

option(USE_MPI "Use MPI" OFF)
if (USE_MPI)
    add_definitions(-DUSE_MPI)
    set(ENABLE_TESTS OFF)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin_mpi" CACHE INTERNAL "Where to place executables")
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/lib_mpi" CACHE INTERNAL "Where to place libraries")
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/lib_mpi" CACHE INTERNAL "Where to place libraries")
    find_package(MPI REQUIRED)
    link_libraries(MPI::MPI_CXX)
endif ()

if (ENABLE_TESTS)
    enable_testing()
    add_subdirectory(Tests)
endif ()


set(TOPLEVEL_INCLUDE_DIR ${PROJECT_SOURCE_DIR})

add_subdirectory(Utilities)
add_subdirectory(Solver)
