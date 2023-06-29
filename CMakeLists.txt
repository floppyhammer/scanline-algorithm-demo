cmake_minimum_required(VERSION 3.16.1)

project(scanline-algorithm-demo)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin")

# Include third_party headers.
include_directories(third_party)

# GLFW library.
include_directories(external/glfw/include) # Include headers.
add_subdirectory(external/glfw)

add_executable(scanline-algorithm-demo
        src/main.cpp
        src/simple_math.cpp)

# Windows.
target_link_libraries(scanline-algorithm-demo glfw)
