add_executable(MathTest DotTest.cpp)
add_executable(SolveTest SolveTest.cpp
        ${TOPLEVEL_INCLUDE_DIR}/solvers.h ${TOPLEVEL_INCLUDE_DIR}/solvers.cpp
        ${TOPLEVEL_INCLUDE_DIR}/csr.h ${TOPLEVEL_INCLUDE_DIR}/csr.cpp
        ${TOPLEVEL_INCLUDE_DIR}/input.h ${TOPLEVEL_INCLUDE_DIR}/input.cpp)

add_dependencies(MathTest MathFunc)
add_dependencies(SolveTest MathFunc)

target_link_libraries(MathTest PRIVATE MathFunc)
target_link_libraries(SolveTest PRIVATE MathFunc)

set_target_properties(MathTest PROPERTIES FOLDER "Tests")
set_target_properties(SolveTest PROPERTIES FOLDER "Tests")


add_test(
        NAME Test1
        COMMAND $<TARGET_FILE:MathTest>
)

add_test(
        NAME Test2
        COMMAND $<TARGET_FILE:SolveTest>
)