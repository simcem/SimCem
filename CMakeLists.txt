cmake_minimum_required(VERSION 3.0)
project(simcem)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

add_subdirectory(src/pybind11)
include_directories(src/stator src/simcem src)


add_library(simcem src/simcem/simcem/simcem.cpp)
target_link_libraries(simcem PUBLIC Eigen3::Eigen)
set_target_properties(simcem PROPERTIES POSITION_INDEPENDENT_CODE True)

pybind11_add_module(core pysrc/simcem/core.cpp)
target_link_libraries(core PUBLIC simcem)
