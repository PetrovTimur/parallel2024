add_executable(MathTests tests.cpp)

add_dependencies(MathTests MathFunc)

target_link_libraries(MathTests PRIVATE MathFunc)

set_target_properties(MathTests PROPERTIES FOLDER "Tests")


add_test(
        NAME MathTest1
        COMMAND $<TARGET_FILE:MathTests>
)