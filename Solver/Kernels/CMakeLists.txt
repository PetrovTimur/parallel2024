add_library(MathFunc SHARED
        mathfunc.cpp
        mathfunc.h)

target_include_directories(MathFunc PUBLIC ${TOPLEVEL_INCLUDE_DIR})


