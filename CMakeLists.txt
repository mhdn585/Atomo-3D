cmake_minimum_required(VERSION 3.10)
project(OpenGLProject)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(glfw3 REQUIRED)
find_package(OpenGL REQUIRED)

# ¡ESTA LÍNEA ES LA IMPORTANTE!
include_directories(${PROJECT_SOURCE_DIR}/external/glad)
include_directories(${PROJECT_SOURCE_DIR}/src)

set(SOURCES
    src/main.cpp
    src/glad.c
)

add_executable(opengl_app ${SOURCES})

target_link_libraries(opengl_app 
    glfw
    ${OPENGL_LIBRARIES}
    ${GLFW_LIBRARIES}
)

target_link_libraries(opengl_app dl)