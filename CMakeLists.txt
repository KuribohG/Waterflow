cmake_minimum_required(VERSION 3.6)
project(FluidSimulation)

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_STANDARD 11)

include(FindOpenGL)
include(FindGLUT)

include_directories(${OPENGL_INCLUDE_DIR} ${GLUT_INCLUDE_DIR})
set(LIBRARIES ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES})

aux_source_directory(FluidSimulation SOURCE_FILES)
add_executable(FluidSimulation ${SOURCE_FILES})
target_link_libraries(FluidSimulation ${LIBRARIES})