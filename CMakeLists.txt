cmake_minimum_required(VERSION 3.15)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_TOOLCHAIN_FILE ./vcpkg/scripts/buildsystems/vcpkg.cmake)
MESSAGE(STATUS "Build type: " ${CMAKE_BUILD_TYPE})
project(simcem)

find_package(PkgConfig)
pkg_check_modules(IPOPT REQUIRED ipopt)
pkg_check_modules(EIGEN3 REQUIRED eigen3)

# if cmake version is higher than 3.20 we can use the new boost configuration
if(POLICY CMP0167)
  cmake_policy(SET CMP0167 NEW) #  Use the new boost configuration
endif()
find_package(Boost 1.70.0 REQUIRED system filesystem program_options unit_test_framework)

find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)
find_package(SUNDIALS REQUIRED COMPONENTS ida)

include_directories(${PROJECT_SOURCE_DIR}/src/stator)

include_directories(src/simcem/ ${EIGEN3_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS} ${IPOPT_INCLUDE_DIRS} ${SUNDIALS_INCLUDE_DIRS} src/stator src/)
link_libraries(${IPOPT_LIBRARIES} ${IPOPT_LINK_FLAGS} ${Boost_LIBRARIES} ${SUNDIALS_LIBRARIES})

add_library(simcem STATIC src/simcem/simcem/simcem.cpp)
set_target_properties(simcem PROPERTIES POSITION_INDEPENDENT_CODE True)

pybind11_add_module(core MODULE pysrc/simcem/core.cpp)
target_link_libraries(core PUBLIC simcem)

pybind11_add_module(kiln MODULE pysrc/simcem/kiln.cpp)
target_link_libraries(kiln PRIVATE simcem SUNDIALS::ida)

#target_include_directories(kiln PUBLIC ${SUNDIALS_INCLUDE_DIRS})

install(TARGETS kiln core DESTINATION ${CMAKE_INSTALL_PREFIX})