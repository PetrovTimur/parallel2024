add_library(MathFunc SHARED
        mathfunc.cpp
        mathfunc.h)

target_include_directories(MathFunc PUBLIC ${PROJECT_SOURCE_DIR})