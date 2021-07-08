cmake_minimum_required(VERSION 3.16.1)
project(scanline-vg-renderer)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin")

include_directories(external/inc)

# Windows
link_directories(external/lib/lib-vc2019)

add_executable(scanline-vg-renderer
        external/inc/glad.c
        src/main.cpp
        src/math/simple.cpp
        src/math/simple.h
        src/nanosvg.h
        )

# Windows
target_link_libraries(scanline-vg-renderer glfw3)
#set(GLFW_INCLUDE_DIR "C:/DevEnv/glfw-3.3.3.bin.WIN64/include")
#set(GLFW_LIBRARY_DIR "C:/DevEnv/glfw-3.3.3.bin.WIN64/lib-vc2019")
#target_include_directories(scanline-vg-renderer PUBLIC ${GLFW_INCLUDE_DIR})
#target_link_libraries(scanline-vg-renderer ${GLFW_LIBRARY_DIR}/glfw3.lib)
#target_link_libraries(scanline-vg-renderer -opengl32)

# Linux
# Do "sudo apt install libglfw3-dev" first
#target_link_libraries(scanline-vg-renderer -lglfw -lGL)
