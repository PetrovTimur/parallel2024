set(MEASURE_SOURCES
    DotMeasure.cpp
    AXpYMeasure.cpp
    SpMVMeasure.cpp
)

foreach (MEASURE_SOURCE ${MEASURE_SOURCES})
    cmake_path(GET MEASURE_SOURCE STEM MEASURE)

    if (${MEASURE} STREQUAL "SpMVMeasure")
        add_executable(${MEASURE} ${MEASURE_SOURCE} ${PROJECT_SOURCE_DIR}/Solver/csr.cpp)
        add_dependencies(${MEASURE} Utilities)
        target_link_libraries(${MEASURE} PUBLIC Utilities)
    else ()
        add_executable(${MEASURE} ${MEASURE_SOURCE})
    endif ()

    add_dependencies(${MEASURE} MathFunc)
    target_link_libraries(${MEASURE} PUBLIC MathFunc)

endforeach ()