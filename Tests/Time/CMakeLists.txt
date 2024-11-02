add_executable(DotMeasure DotMeasure.cpp)
add_executable(AXpYMeasure AXpYMeasure.cpp)

add_dependencies(DotMeasure MathFunc)
add_dependencies(AXpYMeasure MathFunc)

target_link_libraries(DotMeasure PUBLIC MathFunc)
target_link_libraries(AXpYMeasure PUBLIC MathFunc)