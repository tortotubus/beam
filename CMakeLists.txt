cmake_minimum_required(VERSION 3.12.0...4.0.0)
message(STATUS "CMake version: ${CMAKE_VERSION}")

# Require C++17 and disable compiler-specific extensions
set(CMAKE_CXX_STANDARD 17 CACHE STRING "C++ standard to use.")
set(CMAKE_CXX_STANDARD_REQUIRED ON CACHE BOOL
  "Force the use of the chosen C++ standard.")
set(CMAKE_CXX_EXTENSIONS OFF CACHE BOOL "Enable C++ standard extensions.")

project(elff NONE)
# Current version of ELFF, see also `makefile`.
#   mfem_VERSION = (string)
#   ELFF_VERSION = (int)   [automatically derived from mfem_VERSION]
set(${PROJECT_NAME}_VERSION 0.0.1)

# Load project defaults
include("${CMAKE_CURRENT_SOURCE_DIR}/config/defaults.cmake")

# Prohibit in-source build
if (${PROJECT_SOURCE_DIR} STREQUAL ${PROJECT_BINARY_DIR})
  message(FATAL_ERROR
    "ELFF does not support in-source CMake builds at this time.")
endif (${PROJECT_SOURCE_DIR} STREQUAL ${PROJECT_BINARY_DIR})

# Path to CMake modules and utilities.
set(ELFF_CMAKE_PATH ${PROJECT_SOURCE_DIR}/config)
set(CMAKE_MODULE_PATH ${ELFF_CMAKE_PATH}/cmake/modules)

# Load ELFF CMake utilities.
include(ElffCmakeUtilities)

# Enable languages.
enable_language(CXX)

# Suppress warnings about MACOSX_RPATH
set(CMAKE_MACOSX_RPATH OFF CACHE BOOL "")
if (CMAKE_BUILD_TYPE MATCHES "Debug|debug|DEBUG")
  set(ELFF_DEBUG ON)
else()
  set(ELFF_DEBUG OFF)
endif()

# Set 
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

# ELFF_DEBUG
if (CMAKE_BUILD_TYPE MATCHES "Debug|debug|DEBUG")
  set(ELFF_DEBUG ON)
else()
  set(ELFF_DEBUG OFF)
endif()


# BLAS, LAPACK
if (ELFF_USE_LAPACK)
  find_package(BLAS REQUIRED)
  find_package(LAPACK REQUIRED)
endif()

# CoDiPack package
if (MFEM_USE_CODIPACK)
  find_package(CODIPACK REQUIRED)
  # find_package updates CODIPACK_FOUND, CODIPACK_INCLUDE_DIRS, CODIPACK_LIBRARIES
endif()

option (ELFF_USE_MPI "Build elff with MPI" ON)
find_package(MPI REQUIRED)

# Find Eigen3 package
find_package(Eigen3 REQUIRED)

# Find HDF5 package

set(HDF5_PREFER_PARALLEL TRUE)
# set(HDF5_DIR "/usr/lib64/openmpi/lib/cmake/hdf5") 
find_package(HDF5 COMPONENTS C HL REQUIRED)
# find_package(HDF5 COMPONENTS C CXX HL REQUIRED)

# Headers and sources
set(SOURCES "")
set(HEADERS "")
set(ELFF_SOURCE_DIRS general models basilisk io)

# Collect from the specified subdirectories
foreach(DIR IN LISTS ELFF_SOURCE_DIRS)
  add_subdirectory(${DIR})
endforeach()

add_subdirectory(config)

# Variables used when generating _config.hpp, and config.mk
set(ELFF_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(ELFF_INSTALL_DIR ${CMAKE_INSTALL_PREFIX})

# Configure the config header
configure_file(
  "${PROJECT_SOURCE_DIR}/config/cmake/config.hpp.in"
  "${PROJECT_BINARY_DIR}/config/_config.hpp")

# Declaring the library
add_library(elff SHARED ${SOURCES} ${HEADERS})

# --- OpenMP (optional) ---
option(ELFF_USE_OPENMP "Build ELFF with OpenMP" ON)
if (ELFF_USE_OPENMP)
  find_package(OpenMP COMPONENTS CXX REQUIRED)
  target_link_libraries(elff PRIVATE OpenMP::OpenMP_CXX)
  target_compile_definitions(elff PUBLIC ELFF_WITH_OPENMP=1)
endif()

# Add binary dir to include paths to find generated config.hpp
target_include_directories(elff 
  PUBLIC 
    ${CMAKE_CURRENT_BINARY_DIR}/config
    ${CMAKE_CURRENT_SOURCE_DIR}
)

# Link in numeric libs
target_link_libraries(elff 
  PUBLIC
    hdf5::hdf5
    hdf5::hdf5_hl
    MPI::MPI_CXX
    Eigen3::Eigen
    OpenMP::OpenMP_CXX
    # hdf5::hdf5_cpp       # C++ API
)

option(BUILD_TESTS "Build unit & integration tests" ON)
if(BUILD_TESTS)
  # Find GTest
  find_package(GTest REQUIRED)

  enable_testing()
  add_subdirectory(tests)
endif()

# Add Doxygen documentation
option(BUILD_DOC "Build documentation" OFF)
if (BUILD_DOC)
  find_package(Doxygen)
  if (DOXYGEN_FOUND)
      add_custom_target(docs ALL
          ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile
          WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
          COMMENT "Generating API documentation with Doxygen"
          VERBATIM)
  endif()
endif()
