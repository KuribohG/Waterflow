cmake_minimum_required(VERSION 3.5)
project(FluidSimulation)

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_COMPILER g++)
set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} -fopenmp)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()


aux_source_directory(FluidSimulation SOURCE_FILES)
add_executable(FluidSimulation ${SOURCE_FILES})
target_link_libraries(FluidSimulation ${LIBRARIES})
