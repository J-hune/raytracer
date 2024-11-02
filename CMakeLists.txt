cmake_minimum_required(VERSION 3.20)
project(raytracer)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O3")

# Automatic retrieval of source files
file(GLOB_RECURSE SOURCES src/*.cpp)

# Automatic retrieval of header files
file(GLOB_RECURSE HEADERS include/*.h)

# Add executable
add_executable(raytracer ${SOURCES} ${HEADERS})

include_directories(${OPENGL_INCLUDE_DIRS} ${GLUT_INCLUDE_DIRS})
target_link_libraries(raytracer ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES})