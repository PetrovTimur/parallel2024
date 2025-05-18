cmake_minimum_required(VERSION 3.12)
project(parallel2024)

# Set default build type to Release if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif()

find_package(OpenMP)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

# Debug/Release flags
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin" CACHE INTERNAL "Where to place executables")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib" CACHE INTERNAL "Where to place libraries")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib" CACHE INTERNAL "Where to place libraries")
file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/log")      # Where to place logs

#set(TOPLEVEL_INCLUDE_DIR ${PROJECT_SOURCE_DIR})
add_compile_definitions(PROJECT_SOURCE_DIR=\"${PROJECT_SOURCE_DIR}\")
add_compile_definitions(LOG_DIR=\"${PROJECT_BINARY_DIR}/log\")

option(ENABLE_TESTS "Enable automatic testing" ON)
option(DEBUG_MODE "Use debug tools" OFF)
option(USE_MPI "Use MPI" OFF)
option(USE_CUDA "Use CUDA" OFF)
option(MEASURES "Enable measuring" OFF)

if(USE_MPI AND USE_CUDA)
    message(FATAL_ERROR "You cannot enable both USE_MPI and USE_CUDA at the same time.")
endif()

if (DEBUG_MODE)
    add_compile_definitions(DEBUG_MODE)
endif ()

if (USE_MPI)
    add_compile_definitions(USE_MPI)
#    set(ENABLE_TESTS OFF)
    find_package(MPI REQUIRED)
    link_libraries(MPI::MPI_CXX)
endif ()

if (USE_CUDA)
    if(POLICY CMP0104)
        cmake_policy(SET CMP0104 NEW)
    endif()

    if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
        set(CMAKE_CUDA_ARCHITECTURES 89)
    endif()

    add_compile_definitions(USE_CUDA)
    enable_language(CUDA)
#    find_package(CUDAToolkit REQUIRED)
endif ()

if (ENABLE_TESTS)
    enable_testing()
    add_subdirectory(Tests)
endif ()

add_subdirectory(Utilities)
add_subdirectory(Solver)
