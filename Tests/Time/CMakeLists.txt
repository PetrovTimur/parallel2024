set(MEASURES
    DotMeasure
    AXpYMeasure
    SpMVMeasure
    CGSolverMeasure
    CSRAdjacencyMeasure
)

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

endforeach ()