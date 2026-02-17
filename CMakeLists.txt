cmake_minimum_required(VERSION 3.24)
project(raytracer VERSION 0.3.0 LANGUAGES CXX CUDA)

include(FetchContent)

find_package(CUDAToolkit 12 REQUIRED)
find_package(OpenGL REQUIRED)
find_package(PNG REQUIRED)
find_package(ZLIB REQUIRED)
find_package(glfw3 3.4 REQUIRED)

set(FASTGLTF_COMPILE_AS_CPP20 ON CACHE BOOL "" FORCE)
set(FASTGLTF_ENABLE_EXAMPLES OFF CACHE BOOL "" FORCE)
set(FASTGLTF_ENABLE_TESTS OFF CACHE BOOL "" FORCE)
FetchContent_Declare(
    fastgltf
    GIT_REPOSITORY https://github.com/spnda/fastgltf.git
    GIT_TAG v0.9.0
    GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(fastgltf)

FetchContent_Declare(
    stb
    GIT_REPOSITORY https://github.com/nothings/stb.git
    GIT_TAG 31c1ad37456438565541f4919958214b6e762fb4
    GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(stb)

FetchContent_Declare(
    tinyexr
    GIT_REPOSITORY https://github.com/syoyo/tinyexr.git
    GIT_TAG v1.0.13
    GIT_SHALLOW TRUE
    SOURCE_SUBDIR no-cmake-project
)
FetchContent_MakeAvailable(tinyexr)

FetchContent_Declare(
    optix
    GIT_REPOSITORY https://github.com/NVIDIA/optix-sdk.git
    GIT_TAG v9.0.0
    GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(optix)

set(device_ir "${CMAKE_CURRENT_BINARY_DIR}/device.optixir")
set(device_source "${CMAKE_CURRENT_BINARY_DIR}/device_ir.cpp")
add_custom_command(
    OUTPUT "${device_ir}"
    COMMAND "${CMAKE_CUDA_COMPILER}"
        --optix-ir --use_fast_math --std=c++20 --gpu-architecture=compute_75
        "-I${optix_SOURCE_DIR}/include"
        "-I${CUDAToolkit_INCLUDE_DIRS}"
        "-I${CMAKE_CURRENT_SOURCE_DIR}/src"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/device.cu"
        -o "${device_ir}"
    DEPENDS src/device.cu src/gpu_shared.hpp
    VERBATIM
)
add_custom_command(
    OUTPUT "${device_source}"
    COMMAND "${CMAKE_COMMAND}"
        "-DBIN2C=${CUDAToolkit_BIN_DIR}/bin2c"
        "-DINPUT=${device_ir}"
        "-DOUTPUT=${device_source}"
        -P "${CMAKE_CURRENT_SOURCE_DIR}/cmake/embed.cmake"
    DEPENDS "${device_ir}" cmake/embed.cmake
    VERBATIM
)

add_executable(raytracer
    "${device_source}"
    src/display.cpp
    src/image.cpp
    src/main.cpp
    src/photon_grid.cu
    src/renderer.cpp
    src/scene.cpp
    src/stb.cpp
    src/tinyexr.cpp
)
set_source_files_properties(src/tinyexr.cpp PROPERTIES COMPILE_OPTIONS -w)
target_compile_features(raytracer PRIVATE cxx_std_20)
target_compile_options(raytracer PRIVATE
    $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CXX_COMPILER_ID:GNU,Clang>>:-Wall;-Wextra;-Wpedantic>
)
target_include_directories(raytracer PRIVATE
    "${optix_SOURCE_DIR}/include"
    "${stb_SOURCE_DIR}"
    "${tinyexr_SOURCE_DIR}"
)
target_link_libraries(raytracer PRIVATE
    CUDA::cuda_driver
    CUDA::cudart
    OpenGL::GL
    PNG::PNG
    ZLIB::ZLIB
    fastgltf::fastgltf
    glfw
    ${CMAKE_DL_LIBS}
)

set_target_properties(raytracer PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
    CUDA_ARCHITECTURES 75
)
