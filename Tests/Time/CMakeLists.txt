set(MEASURES
    DotMeasure
    AXpYMeasure
    SpMVMeasure
    CGSolverMeasure
)

if (NOT USE_MPI)
    set(MEASURES ${MEASURES} CSRAdjacencyMeasure)
endif ()

foreach (MEASURE ${MEASURES})
#    cmake_path(GET MEASURE_SOURCE STEM MEASURE)
    set(MEASURE_SOURCE ${MEASURE}.cpp)

    if ((${MEASURE} STREQUAL "SpMVMeasure") OR (${MEASURE} STREQUAL "CSRAdjacencyMeasure"))
        add_executable(${MEASURE} ${MEASURE_SOURCE} ${PROJECT_SOURCE_DIR}/Solver/csr.cpp)
        add_dependencies(${MEASURE} IOs)
        target_link_libraries(${MEASURE} PUBLIC IOs)
        if (USE_MPI)
            target_link_libraries(${MEASURE} PUBLIC Coms)
        endif ()
    elseif (${MEASURE} STREQUAL "CGSolverMeasure")
        add_executable(${MEASURE} ${MEASURE_SOURCE} ${PROJECT_SOURCE_DIR}/Solver/csr.cpp ${PROJECT_SOURCE_DIR}/Solver/solvers.cpp)
        add_dependencies(${MEASURE} IOs)
        target_link_libraries(${MEASURE} PUBLIC IOs)
        if (USE_MPI)
            target_link_libraries(${MEASURE} PUBLIC Coms)
        endif ()
    else ()
        add_executable(${MEASURE} ${MEASURE_SOURCE})
    endif ()

    add_dependencies(${MEASURE} MathFunc)
    target_link_libraries(${MEASURE} PUBLIC MathFunc)
    target_include_directories(${MEASURE} PUBLIC ${PROJECT_SOURCE_DIR})

endforeach ()

if (USE_CUDA)
    set(CUDA_MEASURES
            DotMeasureCUDA
            AXpYMeasureCUDA
            SpMVMeasureCUDA
            CGSolverMeasureCUDA
    )

    foreach (MEASURE ${CUDA_MEASURES})
        set(MEASURE_SOURCE ${MEASURE}.cu)

        add_executable(${MEASURE} ${MEASURE_SOURCE} ${PROJECT_SOURCE_DIR}/Solver/csr.cpp ${PROJECT_SOURCE_DIR}/Solver/solvers.cu)
        add_dependencies(${MEASURE} MathFuncCUDA)
        target_link_libraries(${MEASURE} PUBLIC MathFuncCUDA OpenMP::OpenMP_CXX IOs)
        target_include_directories(${MEASURE} PUBLIC ${PROJECT_SOURCE_DIR})
        set_target_properties(${MEASURE} PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
        target_compile_options(${MEASURE} PRIVATE $<$<COMPILE_LANGUAGE:CUDA>: -Xcompiler=-fopenmp>)
        target_compile_features(${MEASURE} PRIVATE cxx_std_14 cuda_std_14)
    endforeach ()
endif ()