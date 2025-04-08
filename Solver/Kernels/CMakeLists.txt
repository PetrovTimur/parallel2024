add_library(MathFunc SHARED
        mathfunc.cpp
        mathfunc.h
        )

target_include_directories(MathFunc PUBLIC ${PROJECT_SOURCE_DIR})


if (USE_CUDA)
        add_library(MathFuncCUDA STATIC
                mathfunc.cu
                mathfunc.cuh)

        target_include_directories(MathFuncCUDA PUBLIC ${PROJECT_SOURCE_DIR} ${CUDAToolkit_INCLUDE_DIRS})
        set_target_properties(MathFuncCUDA PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
endif ()