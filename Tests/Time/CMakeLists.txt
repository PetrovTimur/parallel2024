set(MEASURES
    DotMeasure
    AXpYMeasure
    SpMVMeasure
    CGSolverMeasure
)

foreach (MEASURE ${MEASURES})
#    cmake_path(GET MEASURE_SOURCE STEM MEASURE)
    set(MEASURE_SOURCE ${MEASURE}.cpp)

    if (${MEASURE} STREQUAL "SpMVMeasure")
        add_executable(${MEASURE} ${MEASURE_SOURCE} ${PROJECT_SOURCE_DIR}/Solver/csr.cpp)
        add_dependencies(${MEASURE} IOs)
        target_link_libraries(${MEASURE} PUBLIC IOs Coms)
    elseif (${MEASURE} STREQUAL "CGSolverMeasure")
        add_executable(${MEASURE} ${MEASURE_SOURCE} ${PROJECT_SOURCE_DIR}/Solver/csr.cpp ${PROJECT_SOURCE_DIR}/Solver/solvers.cpp)
        add_dependencies(${MEASURE} IOs)
        target_link_libraries(${MEASURE} PUBLIC IOs Coms)
    else ()
        add_executable(${MEASURE} ${MEASURE_SOURCE})
    endif ()

    add_dependencies(${MEASURE} MathFunc)
    target_link_libraries(${MEASURE} PUBLIC MathFunc)

endforeach ()