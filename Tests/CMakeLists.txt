add_executable(DotCorrectnessTest DotCorrectnessTest.cpp)
add_executable(DotParallelTest DotParallelTest.cpp)
add_executable(AXpYCorrectnessTest AXpYCorrectnessTest.cpp)
add_executable(AXpYParallelTest AXpYParallelTest.cpp)
add_executable(SpMVCorrectnessTest SpMVCorrectnessTest.cpp)

#add_executable(SolveTest SolveTest.cpp
#        ${PROJECT_SOURCE_DIR}/Solver/solvers.h ${TOPLEVEL_INCLUDE_DIR}/solvers.cpp
#        ${TOPLEVEL_INCLUDE_DIR}/csr.h ${TOPLEVEL_INCLUDE_DIR}/csr.cpp
#        ${TOPLEVEL_INCLUDE_DIR}/Utilities/input.h ${TOPLEVEL_INCLUDE_DIR}/Utilities/input.cpp)

add_dependencies(DotCorrectnessTest MathFunc)
add_dependencies(DotParallelTest MathFunc)
add_dependencies(AXpYCorrectnessTest MathFunc)
add_dependencies(AXpYParallelTest MathFunc)
add_dependencies(SpMVCorrectnessTest MathFunc)
#add_dependencies(SolveTest MathFunc)

target_link_libraries(DotCorrectnessTest PRIVATE MathFunc)
target_link_libraries(DotParallelTest PRIVATE MathFunc)
target_link_libraries(AXpYCorrectnessTest PRIVATE MathFunc)
target_link_libraries(AXpYParallelTest PRIVATE MathFunc)
target_link_libraries(SpMVCorrectnessTest PRIVATE MathFunc)
#target_link_libraries(SolveTest PRIVATE MathFunc)

set_target_properties(DotCorrectnessTest PROPERTIES FOLDER "Tests")
set_target_properties(DotParallelTest PROPERTIES FOLDER "Tests")
set_target_properties(AXpYCorrectnessTest PROPERTIES FOLDER "Tests")
set_target_properties(AXpYParallelTest PROPERTIES FOLDER "Tests")
set_target_properties(SpMVCorrectnessTest PROPERTIES FOLDER "Tests")
#set_target_properties(SolveTest PROPERTIES FOLDER "Tests")


add_test(
        NAME DotCorrectnessTest
        COMMAND $<TARGET_FILE:DotCorrectnessTest>
)
add_test(
        NAME DotParallelTest
        COMMAND $<TARGET_FILE:DotParallelTest>
)
add_test(
        NAME AXpYCorrectnessTest
        COMMAND $<TARGET_FILE:AXpYCorrectnessTest>
)
add_test(
        NAME AXpYParallelTest
        COMMAND $<TARGET_FILE:AXpYParallelTest>
)
add_test(
        NAME SpMVCorrectnessTest.cpp
        COMMAND $<TARGET_FILE:SpMVCorrectnessTest>
)


#add_test(
#        NAME Test2
#        COMMAND $<TARGET_FILE:SolveTest>
#)