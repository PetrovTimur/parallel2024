set(TESTS
        DotCorrectnessTest
        DotParallelTest
        AXpYCorrectnessTest
        AXpYParallelTest
        SpMVCorrectnessTest
        CSRTransposeTest
        CSRBuildAdjacencyMatrixTest
        CSRIncidenceTest
        ReadDataTest
)

if (NOT USE_MPI)
    foreach (TEST ${TESTS})
        #    cmake_path(GET MEASURE_SOURCE STEM MEASURE)
        set(TEST_SOURCE ${TEST}.cpp)

        if (${TEST} MATCHES ".*CSR.*" OR ${TEST} MATCHES ".*Data.*")
            add_executable(${TEST} ${TEST_SOURCE} ${PROJECT_SOURCE_DIR}/Solver/csr.cpp ${PROJECT_SOURCE_DIR}/Utilities/input.cpp)
            target_include_directories(${TEST} PRIVATE ${PROJECT_SOURCE_DIR}/Solver)
        else ()
            add_executable(${TEST} ${TEST_SOURCE})
        endif ()

        # add_executable(${TEST} ${TEST_SOURCE})
        add_dependencies(${TEST} MathFunc)
        target_link_libraries(${TEST} PUBLIC MathFunc)
        set_target_properties(${TEST} PROPERTIES FOLDER "Tests")
        add_test(
                NAME ${TEST}
                COMMAND $<TARGET_FILE:${TEST}>
        )
        if (${TEST} MATCHES ".*Parallel.*")
            set_property(TEST ${TEST} PROPERTY ENVIRONMENT "OMP_NUM_THREADS=4")
        endif ()
    endforeach ()
endif ()

if (USE_CUDA)
    message("CUDA TESTING")
    set(CUDA_TESTS
            AXpYCUDATest
            ReduceCUDATest
            DotCUDATest
            DeviceQueryTest
    )

    foreach (TEST ${CUDA_TESTS})
        set(TEST_SOURCE ${TEST}.cu)
        add_executable(${TEST} ${TEST_SOURCE})
        add_dependencies(${TEST} MathFuncCUDA)
        target_link_libraries(${TEST} PUBLIC MathFuncCUDA OpenMP::OpenMP_CXX)
        set_target_properties(${TEST} PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
        target_compile_options(${TEST} PRIVATE $<$<COMPILE_LANGUAGE:CUDA>: -Xcompiler=-fopenmp>)
        target_compile_features(${TEST} PRIVATE cxx_std_17 cuda_std_17)
    endforeach ()
endif ()

if (MEASURES)
    add_subdirectory(Time)
endif ()
