set(TESTS
        DotCorrectnessTest
        DotParallelTest
        AXpYCorrectnessTest
        AXpYParallelTest
        SpMVCorrectnessTest
        CSRTransposeTest
        CSRBuildAdjacencyMatrixTest
)

if (NOT USE_MPI)
    foreach (TEST ${TESTS})
        #    cmake_path(GET MEASURE_SOURCE STEM MEASURE)
        set(TEST_SOURCE ${TEST}.cpp)

        if (${TEST} MATCHES ".*CSR.*")
            add_executable(${TEST} ${TEST_SOURCE} ${PROJECT_SOURCE_DIR}/Solver/csr.cpp)
            target_include_directories(${TEST} PRIVATE ${PROJECT_SOURCE_DIR}/Solver)
        else ()
            add_executable(${TEST} ${TEST_SOURCE})
        endif ()

#        add_executable(${TEST} ${TEST_SOURCE})
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


#add_executable(SolveTest SolveTest.cpp
#        ${PROJECT_SOURCE_DIR}/Solver/solvers.h ${TOPLEVEL_INCLUDE_DIR}/solvers.cpp
#        ${TOPLEVEL_INCLUDE_DIR}/csr.h ${TOPLEVEL_INCLUDE_DIR}/csr.cpp
#        ${TOPLEVEL_INCLUDE_DIR}/Utilities/input.h ${TOPLEVEL_INCLUDE_DIR}/Utilities/input.cpp)

#add_dependencies(SolveTest MathFunc)
#target_link_libraries(SolveTest PRIVATE MathFunc)
#set_target_properties(SolveTest PROPERTIES FOLDER "Tests")
#add_test(
#        NAME Test2
#        COMMAND $<TARGET_FILE:SolveTest>
#)

add_subdirectory(Time)