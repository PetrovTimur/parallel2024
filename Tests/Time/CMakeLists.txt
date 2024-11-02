add_executable(DotMeasure DotMeasure.cpp)


add_dependencies(DotMeasure MathFunc)

target_link_libraries(DotMeasure PUBLIC MathFunc)