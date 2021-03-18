cmake_minimum_required(VERSION 3.0)
  
project(simcem)

add_subdirectory(src/pybind11)

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
find_package(IPOPT 3.11.9 REQUIRED)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(Boost REQUIRED )
find_package(SUNDIALS REQUIRED )

message(STATUS "CMakeLists: Checking ${IPOPT_INCLUDE_DIRS}")

include_directories(src/stator src/simcem src)
add_library(simcem src/simcem/simcem/simcem.cpp)

target_include_directories(simcem PUBLIC src/simcem/ ${EIGEN3_INCLUDE_DIR} ${IPOPT_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})
target_link_libraries(simcem PUBLIC PRIVATE Eigen3::Eigen ${IPOPT_LIBRARIES} ${IPOPT_LINK_FLAGS})
set_target_properties(simcem PROPERTIES POSITION_INDEPENDENT_CODE True)

pybind11_add_module(core pysrc/simcem/core.cpp)
target_link_libraries(core PRIVATE simcem)

pybind11_add_module(kiln pysrc/simcem/kiln.cpp)
target_link_libraries(kiln PRIVATE simcem ${SUNDIALS_LIBRARIES})
target_include_directories(kiln PUBLIC ${SUNDIALS_INCLUDE_DIRS})
