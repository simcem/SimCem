cmake_minimum_required(VERSION 3.15)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# If we're building on nix, then skip using vcpkg
if(DEFINED ENV{NIX_BUILD_TOP})
  message(STATUS "Detected Nix build environment, skipping vcpkg dependency management.")
else()
  set(CMAKE_TOOLCHAIN_FILE ./submodule/vcpkg/scripts/buildsystems/vcpkg.cmake)
endif()
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


set(COINHSL_EXTRACT_NAME "coinhsl-archive-2024.05.15")
set(COINHSL_ARCHIVE_FILE "${CMAKE_SOURCE_DIR}/${COINHSL_EXTRACT_NAME}.tar")
set(COINHSL_EXTRACT_DIR "${CMAKE_SOURCE_DIR}/submodule/ThirdParty-HSL")

# If the COINHSL_DOWNLOAD_URL environmental variable is set, use it to download the archive
if(DEFINED ENV{COINHSL_DOWNLOAD_URL})
  set(COINHSL_DOWNLOAD_URL $ENV{COINHSL_DOWNLOAD_URL})
  message(STATUS "Downloading CoinHSL archive from ${COINHSL_DOWNLOAD_URL}")
  file(DOWNLOAD "${COINHSL_DOWNLOAD_URL}" "${COINHSL_ARCHIVE_FILE}" SHOW_PROGRESS STATUS download_status)
  list(GET download_status 0 download_result)
  if(NOT download_result EQUAL 0)
    message(FATAL_ERROR "Failed to download CoinHSL archive from ${COINHSL_DOWNLOAD_URL}")
  endif()
endif()

if(NOT EXISTS "${COINHSL_ARCHIVE_FILE}")
  message(WARNING "Skipping build of libcoinhsl.so (which is required). Download ${COINHSL_ARCHIVE_FILE}.tar from https://licences.stfc.ac.uk/product/coin-hsl-archive and place in the root of the source tree if you want to build it.")
else()
  # Check if already extracted?
  if(NOT EXISTS "${COINHSL_EXTRACT_DIR}/coinhsl")
    # Extract the archive
    message(STATUS "Extracting CoinHSL archive to ${COINHSL_EXTRACT_DIR}")
    file(ARCHIVE_EXTRACT INPUT "${COINHSL_ARCHIVE_FILE}" DESTINATION "${COINHSL_EXTRACT_DIR}")
    
    # Rename extracted directory to 'coinhsl'
    if(EXISTS "${COINHSL_EXTRACT_DIR}/${COINHSL_EXTRACT_NAME}")
      message(STATUS "Renaming ${COINHSL_EXTRACT_NAME} to coinhsl")
      file(RENAME "${COINHSL_EXTRACT_DIR}/${COINHSL_EXTRACT_NAME}" "${COINHSL_EXTRACT_DIR}/coinhsl")
    elseif(NOT EXISTS "${COINHSL_EXTRACT_DIR}/coinhsl")
      message(FATAL_ERROR "Could not find extracted directory. Expected either '${COINHSL_EXTRACT_DIR}/${COINHSL_EXTRACT_NAME}' or '${COINHSL_EXTRACT_DIR}/coinhsl'")
    endif()

    message(STATUS "CoinHSL archive extracted and ready")
  endif()

  include(ExternalProject)
  ExternalProject_Add(
    libcoinhsl
    SOURCE_DIR ${PROJECT_SOURCE_DIR}/submodule/ThirdParty-HSL
    CONFIGURE_COMMAND 
      ${PROJECT_SOURCE_DIR}/submodule/ThirdParty-HSL/configure  --prefix=${CMAKE_CURRENT_BINARY_DIR} --libdir=${CMAKE_CURRENT_BINARY_DIR} 
    PREFIX ${PROJECT_SOURCE_DIR}/submodule/ThirdParty-HSL/
    BUILD_COMMAND make
    BUILD_IN_SOURCE 1
  )

  # Create symbolic link from libcoinhsl.so to libhsl.so in the build directory
  add_custom_command(TARGET libcoinhsl POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${CMAKE_CURRENT_BINARY_DIR}/libcoinhsl.so
      ${CMAKE_CURRENT_BINARY_DIR}/libhsl.so
  )

  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libhsl.so DESTINATION ${CMAKE_INSTALL_PREFIX})
endif()


add_library(simcem STATIC src/simcem/simcem/simcem.cpp)
set_target_properties(simcem PROPERTIES POSITION_INDEPENDENT_CODE True)
target_link_libraries(simcem PUBLIC Boost::filesystem Boost::program_options PkgConfig::IPOPT PkgConfig::EIGEN3)
target_include_directories(simcem PUBLIC ${PROJECT_SOURCE_DIR}/src/stator ${PROJECT_SOURCE_DIR}/src/simcem/ ${PROJECT_SOURCE_DIR}/src/)
install(FILES pysrc/simcem/free_database.xml DESTINATION ${CMAKE_INSTALL_PREFIX})

pybind11_add_module(core MODULE pysrc/simcem/core.cpp)
target_link_libraries(core PUBLIC simcem) 

pybind11_add_module(kiln_core MODULE pysrc/simcem/kiln_core.cpp)
target_link_libraries(kiln_core PRIVATE simcem SUNDIALS::ida)

install(TARGETS kiln_core core DESTINATION ${CMAKE_INSTALL_PREFIX})

enable_testing()
add_executable(simcem_example tests/simcem_example.cpp)
target_link_libraries(simcem_example PRIVATE simcem)
add_test(NAME simcem_example_test COMMAND simcem_example)