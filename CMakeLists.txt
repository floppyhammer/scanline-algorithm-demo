cmake_minimum_required(VERSION 3.16.1)
project(scanline-vg-renderer)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin")

# Include external and third_party headers.
include_directories(external third_party)

# GLFW library.
include_directories(external/glfw/include) # Include headers.
add_subdirectory(external/glfw)

add_executable(scanline-vg-renderer
        external/glad.c
        src/main.cpp
        src/math/simple.cpp
        src/math/simple.h
        third_party/nanosvg.h)

# Windows
target_link_libraries(scanline-vg-renderer glfw)
