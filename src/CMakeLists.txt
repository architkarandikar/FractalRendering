cmake_minimum_required(VERSION 2.8)
project(Project1)

# Check compiler
# UNIX
if(UNIX)
  set(GCC_CXX_FLAGS "-std=gnu++11 -m64 -O3 -funroll-loops")
  set(GCC_CXX_FLAGS "${GCC_CXX_FLAGS} -fopenmp")
endif(UNIX)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_CXX_FLAGS}")

# Add modules
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules/")

# Find dependencies
find_package(OpenGL REQUIRED)
find_package(Freetype REQUIRED)
find_package(GLEW REQUIRED)
find_package(CMU462 REQUIRED)

# Install settings
set(CMAKE_INSTALL_PREFIX "${PROJECT_SOURCE_DIR}/../")

# Ubuntu
find_package(PkgConfig REQUIRED)
pkg_search_module(GLFW REQUIRED glfw3)

# Set include directories
include_directories(
    ${FREETYPE_INCLUDE_DIRS}
    ${CMU462_INCLUDE_DIRS}
    ${GLEW_INCLUDE_DIRS}
    ${GLFW_INCLUDE_DIRS}
)

# Set link directories
link_directories(
    ${FREETYPE_LIBRARY_DIRS}
    ${CMU462_LIBRARY_DIRS}
    ${GLEW_LIBRARY_DIRS}
    ${GLFW_LIBRARY_DIRS}
)

# Set drawsvg source
set(CMU462_DRAWSVG_SOURCE
    svg.cpp
    png.cpp
    texture.cpp
    viewport.cpp
    triangulation.cpp
#    hardware_renderer.cpp
    software_renderer.cpp
    drawsvg.cpp
    main.cpp
)

# Set drawsvg header
set(CMU462_DRAWSVG_HEADER
    svg.h
    png.h
    texture.h
    viewport.h
    triangulation.h
#    hardware_renderer.h
    software_renderer.h
    drawsvg.h
)

# Import drawsvg reference
include(reference/reference.cmake)

# drawsvg executable
add_executable( draw_svg
    ${CMU462_DRAWSVG_SOURCE}
    ${CMU462_DRAWSVG_HEADER}
)

# Link drawsvg executable (static link reference solution)
target_link_libraries( draw_svg drawsvg_ref
    ${CMU462_LIBRARIES}
    ${FREETYPE_LIBRARIES}
    ${OPENGL_LIBRARIES}
    ${GLEW_LIBRARIES}
    ${GLFW_LIBRARIES}
)

install(TARGETS draw_svg DESTINATION .)
