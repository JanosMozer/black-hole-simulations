cmake_minimum_required(VERSION 3.21)
project(BlackHoleSim)

#find_package(OpenMP)
#if(OpenMP_FOUND)
#    message(STATUS "OpenMP found, enabling parallelization for CPU fallback.")
#endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Dependencies
find_package(GLEW REQUIRED)
find_package(glfw3 CONFIG REQUIRED)
find_package(glm CONFIG REQUIRED)

# Common dependencies
set(DEPS glfw GLEW::GLEW glm::glm)

# 2D Lensing executable
add_executable(BlackHole2D 2D_lensing.cpp)
target_link_libraries(BlackHole2D PRIVATE ${DEPS})
target_include_directories(BlackHole2D PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_definitions(BlackHole2D PRIVATE GLFW_COCOA_GL_LOADER)
target_compile_options(BlackHole2D PRIVATE -Wno-deprecated-declarations)

# 3D Black Hole executable
add_executable(BlackHole3D black_hole.cpp)
target_link_libraries(BlackHole3D PRIVATE ${DEPS})
#if(OpenMP_FOUND)
#    target_link_libraries(BlackHole3D PRIVATE OpenMP::OpenMP_CXX)
#endif()
target_include_directories(BlackHole3D PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# Simple Lensing executable (OpenGL 4.1 compatible)
# add_executable(SimpleLensing simple_lensing.cpp)
# target_link_libraries(SimpleLensing PRIVATE ${DEPS})
# target_include_directories(SimpleLensing PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# Shader files (copy to output dir)
file(GLOB SHADERS
    "${CMAKE_CURRENT_SOURCE_DIR}/*.vert"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.frag"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.comp"
)

foreach(shader ${SHADERS})
    add_custom_command(TARGET BlackHole3D POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
        ${shader}
        $<TARGET_FILE_DIR:BlackHole3D>
    )
endforeach()

