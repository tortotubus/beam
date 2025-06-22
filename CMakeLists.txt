cmake_minimum_required(VERSION 3.12.0...4.0.0)
message(STATUS "CMake version: ${CMAKE_VERSION}")

# Require C++17 and disable compiler-specific extensions
set(CMAKE_CXX_STANDARD 17 CACHE STRING "C++ standard to use.")
set(CMAKE_CXX_STANDARD_REQUIRED ON CACHE BOOL
  "Force the use of the chosen C++ standard.")
set(CMAKE_CXX_EXTENSIONS OFF CACHE BOOL "Enable C++ standard extensions.")

project(beam NONE)
# Current version of BEAM, see also `makefile`.
#   mfem_VERSION = (string)
#   BEAM_VERSION = (int)   [automatically derived from mfem_VERSION]
set(${PROJECT_NAME}_VERSION 0.0.1)

# Load project defaults
include("${CMAKE_CURRENT_SOURCE_DIR}/config/defaults.cmake")

# Prohibit in-source build
if (${PROJECT_SOURCE_DIR} STREQUAL ${PROJECT_BINARY_DIR})
  message(FATAL_ERROR
    "BEAM does not support in-source CMake builds at this time.")
endif (${PROJECT_SOURCE_DIR} STREQUAL ${PROJECT_BINARY_DIR})

# Path to CMake modules and utilities.
set(BEAM_CMAKE_PATH ${PROJECT_SOURCE_DIR}/config)
set(CMAKE_MODULE_PATH ${BEAM_CMAKE_PATH}/cmake/modules)

# Load BEAM CMake utilities.
include(BeamCmakeUtilities)

# Enable languages.
enable_language(CXX)

# Suppress warnings about MACOSX_RPATH
set(CMAKE_MACOSX_RPATH OFF CACHE BOOL "")
if (CMAKE_BUILD_TYPE MATCHES "Debug|debug|DEBUG")
  set(BEAM_DEBUG ON)
else()
  set(BEAM_DEBUG OFF)
endif()

# Set 
if (BEAM_PRECISION MATCHES "^(double|Double|DOUBLE)$")
  set(BEAM_USE_SINGLE OFF)
  set(BEAM_USE_DOUBLE ON)
elseif (BEAM_PRECISION MATCHES "^(single|Single|SINGLE)$")
  set(BEAM_USE_SINGLE ON)
  set(BEAM_USE_DOUBLE OFF)
else()
  message(FATAL_ERROR " *** Invalid floating-point precision: "
    "BEAM_PRECISION = ${BEAM_PRECISION}")
endif()

# BEAM_DEBUG
if (CMAKE_BUILD_TYPE MATCHES "Debug|debug|DEBUG")
  set(BEAM_DEBUG ON)
else()
  set(BEAM_DEBUG OFF)
endif()

# BLAS, LAPACK
if (BEAM_USE_LAPACK)
  find_package(BLAS REQUIRED)
  find_package(LAPACK REQUIRED)
endif()

# CoDiPack package
if (MFEM_USE_CODIPACK)
  find_package(CODIPACK REQUIRED)
  # find_package updates CODIPACK_FOUND, CODIPACK_INCLUDE_DIRS, CODIPACK_LIBRARIES
endif()

# Find Eigen3 package
find_package(Eigen3 REQUIRED)


# Headers and sources
set(SOURCES "")
set(HEADERS "")
set(BEAM_SOURCE_DIRS general models basilisk)

# Collect from the specified subdirectories
foreach(DIR IN LISTS BEAM_SOURCE_DIRS)
  add_subdirectory(${DIR})
endforeach()

add_subdirectory(config)

# Variables used when generating _config.hpp, and config.mk
set(BEAM_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(BEAM_INSTALL_DIR ${CMAKE_INSTALL_PREFIX})

# Configure the config header
configure_file(
  "${PROJECT_SOURCE_DIR}/config/cmake/config.hpp.in"
  "${PROJECT_BINARY_DIR}/config/_config.hpp")

# Declaring the library
add_library(beam SHARED ${SOURCES} ${HEADERS})

# Add binary dir to include paths to find generated config.hpp
target_include_directories(beam 
  PUBLIC 
    ${CMAKE_CURRENT_BINARY_DIR}/config
    ${CMAKE_CURRENT_SOURCE_DIR}
)

# Link in numeric libs
target_link_libraries(beam 
  PUBLIC
    Eigen3::Eigen
)

option(BUILD_TESTS "Build unit & integration tests" OFF)
if(BUILD_TESTS)
  # Find GTest
  find_package(GTest REQUIRED)

  enable_testing()
  add_subdirectory(tests)
endif()

# Add Doxygen documentation
option(BUILD_DOC "Build documentation" ON)
if (BUILD_DOC)
  find_package(Doxygen)
  if (DOXYGEN_FOUND)
      add_custom_target(doc ALL
          ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile
          WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
          COMMENT "Generating API documentation with Doxygen"
          VERBATIM)
  endif()
endif()
