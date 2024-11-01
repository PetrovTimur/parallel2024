add_executable(DotCorrectnessTest DotCorrectnessTest.cpp)
add_executable(DotParallelTest DotParallelTest.cpp)
#add_executable(SolveTest SolveTest.cpp
#        ${PROJECT_SOURCE_DIR}/Solver/solvers.h ${TOPLEVEL_INCLUDE_DIR}/solvers.cpp
#        ${TOPLEVEL_INCLUDE_DIR}/csr.h ${TOPLEVEL_INCLUDE_DIR}/csr.cpp
#        ${TOPLEVEL_INCLUDE_DIR}/Utilities/input.h ${TOPLEVEL_INCLUDE_DIR}/Utilities/input.cpp)

add_dependencies(DotCorrectnessTest MathFunc)
add_dependencies(DotParallelTest MathFunc)
#add_dependencies(SolveTest MathFunc)

target_link_libraries(DotCorrectnessTest PRIVATE MathFunc)
target_link_libraries(DotParallelTest PRIVATE MathFunc)
#target_link_libraries(SolveTest PRIVATE MathFunc)

set_target_properties(DotCorrectnessTest PROPERTIES FOLDER "Tests")
set_target_properties(DotParallelTest PROPERTIES FOLDER "Tests")
#set_target_properties(SolveTest PROPERTIES FOLDER "Tests")


add_test(
        NAME DotCorrectnessTest
        COMMAND $<TARGET_FILE:DotCorrectnessTest>
)
add_test(
        NAME DotParallelTest
        COMMAND $<TARGET_FILE:DotParallelTest>
)


#add_test(
#        NAME Test2
#        COMMAND $<TARGET_FILE:SolveTest>
#)