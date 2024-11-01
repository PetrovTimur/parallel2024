add_executable(DotTest DotCorrectnessTest.cpp)
#add_executable(SolveTest SolveTest.cpp
#        ${PROJECT_SOURCE_DIR}/Solver/solvers.h ${TOPLEVEL_INCLUDE_DIR}/solvers.cpp
#        ${TOPLEVEL_INCLUDE_DIR}/csr.h ${TOPLEVEL_INCLUDE_DIR}/csr.cpp
#        ${TOPLEVEL_INCLUDE_DIR}/Utilities/input.h ${TOPLEVEL_INCLUDE_DIR}/Utilities/input.cpp)

add_dependencies(DotTest MathFunc)
#add_dependencies(SolveTest MathFunc)

target_link_libraries(DotTest PRIVATE MathFunc)
#target_link_libraries(SolveTest PRIVATE MathFunc)

set_target_properties(DotTest PROPERTIES FOLDER "Tests")
#set_target_properties(SolveTest PROPERTIES FOLDER "Tests")


add_test(
        NAME DotCorrectnessTest
        COMMAND $<TARGET_FILE:DotTest>
)

#add_test(
#        NAME Test2
#        COMMAND $<TARGET_FILE:SolveTest>
#)