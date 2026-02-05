cmake_minimum_required(VERSION 3.24)
project(raytracer VERSION 0.2.0 LANGUAGES CXX)

include(FetchContent)

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

add_executable(raytracer
    src/main.cpp
    src/scene.cpp
)
target_compile_features(raytracer PRIVATE cxx_std_20)
target_compile_options(raytracer PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang>:-Wall;-Wextra;-Wpedantic>
)
target_link_libraries(raytracer PRIVATE fastgltf::fastgltf)

set_target_properties(raytracer PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)
