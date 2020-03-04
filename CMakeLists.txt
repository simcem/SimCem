cmake_minimum_required(VERSION 3.0)
  
project(simcem)

add_subdirectory(src/pybind11)

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
find_package(IPOPT 3.11.9 REQUIRED)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

include_directories(src/stator src/simcem src)
add_library(simcem src/simcem/simcem/simcem.cpp)

target_include_directories(simcem PUBLIC ${EIGEN3_INCLUDE_DIR} ${IPOPT_INCLUDE_DIRS})
target_link_libraries(simcem PUBLIC PRIVATE Eigen3::Eigen ${IPOPT_LIBRARIES} ${IPOPT_LINK_FLAGS} )
set_target_properties(simcem PROPERTIES POSITION_INDEPENDENT_CODE True)

pybind11_add_module(core pysrc/simcem/core.cpp)
target_link_libraries(core PRIVATE simcem)
