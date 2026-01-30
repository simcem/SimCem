cmake_minimum_required(VERSION 3.15)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_TOOLCHAIN_FILE ./submodule/vcpkg/scripts/buildsystems/vcpkg.cmake)
MESSAGE(STATUS "Build type: " ${CMAKE_BUILD_TYPE})
project(simcem)

add_subdirectory(submodule/pybind11)

find_package(PkgConfig)
pkg_check_modules(IPOPT REQUIRED IMPORTED_TARGET GLOBAL ipopt)
pkg_check_modules(EIGEN3 REQUIRED IMPORTED_TARGET GLOBAL eigen3)

# if cmake version is higher than 3.20 we can use the new boost configuration
if(POLICY CMP0167)
  cmake_policy(SET CMP0167 NEW) #  Use the new boost configuration
endif()
find_package(Boost 1.70.0 REQUIRED system filesystem program_options unit_test_framework math)
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(SUNDIALS REQUIRED COMPONENTS ida)

find_package(LAPACK REQUIRED)

include(ExternalProject)

# Custom target to download and extract HSL archive
set(COINHSL_ARCHIVE_URL "https://cloud.bansci.com/public.php/dav/files/iRe55SEGaZB5fZj")
set(COINHSL_ARCHIVE_FILE "${PROJECT_SOURCE_DIR}/coinhsl-archive-2024.05.15.tar")
set(COINHSL_EXTRACT_DIR "${PROJECT_SOURCE_DIR}/submodule/ThirdParty-HSL")
set(COINHSL_EXTRACT_NAME "coinhsl-archive-2024.05.15")

string(REPLACE ";" " " AUTOTOOLS_LAPACK_LIBS "${LAPACK_LIBRARIES}")
message(STATUS "LAPACK found: ${AUTOTOOLS_LAPACK_LIBS}")
add_custom_target(
  download-extract-coinhsl
  COMMENT "Downloading and extracting CoinHSL archive"
  COMMAND ${CMAKE_COMMAND} -P ${PROJECT_SOURCE_DIR}/cmake/DownloadCoinhsl.cmake 
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
  BYPRODUCTS ${COINHSL_EXTRACT_DIR}/coinhsl
)

ExternalProject_Add(
  libcoinhsl
  SOURCE_DIR ${PROJECT_SOURCE_DIR}/submodule/ThirdParty-HSL
  CONFIGURE_COMMAND 
    ${PROJECT_SOURCE_DIR}/submodule/ThirdParty-HSL/configure  --prefix=${CMAKE_CURRENT_BINARY_DIR} --libdir=${CMAKE_CURRENT_BINARY_DIR} 
  PREFIX ${PROJECT_SOURCE_DIR}/submodule/ThirdParty-HSL/
  BUILD_COMMAND make
  BUILD_IN_SOURCE 1
  DEPENDS download-extract-coinhsl
)

# Create symbolic link from libcoinhsl.so to libhsl.so in the build directory
add_custom_command(TARGET libcoinhsl POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E create_symlink
    ${CMAKE_CURRENT_BINARY_DIR}/libcoinhsl.so
    ${CMAKE_CURRENT_BINARY_DIR}/libhsl.so
)

install(FILES 
  ${CMAKE_CURRENT_BINARY_DIR}/libhsl.so
  pysrc/simcem/free_database.xml 
  DESTINATION ${CMAKE_INSTALL_PREFIX}
)

#include_directories(${PROJECT_SOURCE_DIR}/src/stator)
#include_directories(src/simcem/ ${EIGEN3_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS} ${IPOPT_INCLUDE_DIRS} ${SUNDIALS_INCLUDE_DIRS} src/stator src/)

add_library(simcem STATIC src/simcem/simcem/simcem.cpp)
set_target_properties(simcem PROPERTIES POSITION_INDEPENDENT_CODE True)
target_link_libraries(simcem PUBLIC Boost::filesystem Boost::program_options PkgConfig::IPOPT PkgConfig::EIGEN3)
target_include_directories(simcem PUBLIC ${PROJECT_SOURCE_DIR}/src/stator ${PROJECT_SOURCE_DIR}/src/simcem/ ${PROJECT_SOURCE_DIR}/src/)

pybind11_add_module(core MODULE pysrc/simcem/core.cpp)
target_link_libraries(core PUBLIC simcem) 

pybind11_add_module(kiln MODULE pysrc/simcem/kiln.cpp)
target_link_libraries(kiln PRIVATE simcem SUNDIALS::ida)

install(TARGETS kiln core DESTINATION ${CMAKE_INSTALL_PREFIX})

enable_testing()
add_executable(simcem_example tests/simcem_example.cpp)
target_link_libraries(simcem_example PRIVATE simcem)
add_test(NAME simcem_example_test COMMAND simcem_example)