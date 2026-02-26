#
# CMake Settings
# 
cmake_minimum_required(VERSION 3.12.0...4.0.0)
message(STATUS "CMake version: ${CMAKE_VERSION}")

#
# Project Names
#
project(elff NONE)
set(${PROJECT_NAME}_VERSION 0.0.1)

#
# Enable languages
#
enable_language(C)
enable_language(CXX)

#
# C/C++ Language Settings
#
set(CMAKE_CXX_STANDARD 17 CACHE STRING "C++ standard to use.")
set(CMAKE_CXX_STANDARD_REQUIRED ON CACHE BOOL
  "Force the use of the chosen C++ standard.")
set(CMAKE_CXX_EXTENSIONS OFF CACHE BOOL "Enable C++ standard extensions.")

#
# Default Options
#
include("${CMAKE_CURRENT_SOURCE_DIR}/config/defaults.cmake")

#
# Set CMAKE_MODULE_PATH
#
set(ELFF_CMAKE_PATH ${PROJECT_SOURCE_DIR}/config)
set(CMAKE_MODULE_PATH ${ELFF_CMAKE_PATH}/cmake/modules)

# 
# Load ELFF CMake utilities.
#
include(ElffCmakeUtilities)

# 
# Floating point precision settings
#
if (ELFF_PRECISION MATCHES "^(double|Double|DOUBLE)$")
  set(ELFF_USE_SINGLE OFF)
  set(ELFF_USE_DOUBLE ON)
elseif (ELFF_PRECISION MATCHES "^(single|Single|SINGLE)$")
  set(ELFF_USE_SINGLE ON)
  set(ELFF_USE_DOUBLE OFF)
else()
  message(FATAL_ERROR " *** Invalid floating-point precision: "
    "ELFF_PRECISION = ${ELFF_PRECISION}")
endif()

#
# Collect project sources
# 
set(SOURCES "")
set(HEADERS "")
set(ELFF_SOURCE_DIRS c config general io models mpi)

foreach(DIR IN LISTS ELFF_SOURCE_DIRS)
  add_subdirectory(${DIR})
endforeach()


# 
# Configure config/_config.hpp which is included in config/config.hpp
#
configure_file(
  "${PROJECT_SOURCE_DIR}/config/cmake/config.hpp.in"
  "${PROJECT_BINARY_DIR}/config/_config.hpp")


#
# Declare shared (.so) library
#
add_library(elff SHARED ${SOURCES} ${HEADERS})

#
# OpenMP
#
if (ELFF_USE_OPENMP)
  find_package(OpenMP COMPONENTS CXX REQUIRED)
  target_link_libraries(elff PUBLIC OpenMP::OpenMP_CXX)
  target_compile_definitions(elff PUBLIC ELFF_WITH_OPENMP=1)
endif()

#
# MPI
#
if (ELFF_USE_MPI)
  find_package(MPI COMPONENTS C CXX REQUIRED)
  target_link_libraries(elff PUBLIC MPI::MPI_CXX)
endif()

#
# Eigen3
#
if (ELFF_USE_EIGEN3)
  find_package(Eigen3 REQUIRED)
  target_link_libraries(elff PUBLIC Eigen3::Eigen)
endif()

#
# HDF5
#
if (ELFF_USE_HDF5)
  if (ELFF_USE_MPI)
    set(HDF5_PREFER_PARALLEL TRUE)
  endif()
  find_package(HDF5 COMPONENTS C HL REQUIRED)
  target_link_libraries(elff PUBLIC
    hdf5::hdf5
    hdf5::hdf5_hl
  )
endif()

#
# Add ${BINARY_DIR}/config to include path so that configured file ${BINARY_DIR}/config/_config.hpp is found
#
target_include_directories(elff 
  PUBLIC 
    ${CMAKE_CURRENT_BINARY_DIR}/config
    ${CMAKE_CURRENT_SOURCE_DIR}
)

#
# Basilisk
#
if (ELFF_USE_BASILISK)
  add_subdirectory(basilisk)
endif()

#
# Tests
#
if(ELFF_BUILD_TESTS)
  # Find GTest
  find_package(GTest REQUIRED)

  enable_testing()
  add_subdirectory(tests)
endif()

#
# Documentation
#
if (ELFF_BUILD_DOC)
  find_package(Doxygen)
  if (DOXYGEN_FOUND)
      add_custom_target(docs ALL
          ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile
          WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
          COMMENT "Generating API documentation with Doxygen"
          VERBATIM)
  endif()
endif()
