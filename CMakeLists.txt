#
# CMake Settings
# 
cmake_minimum_required(VERSION 3.12.0...4.0.0)
message(STATUS "CMake version: ${CMAKE_VERSION}")

#
# Project Names
#
project(ELFF
  VERSION 0.0.1
)
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

include(GNUInstallDirs)

if(DEFINED CMAKE_MPI_INSTALL_LIBDIR AND CMAKE_MPI_INSTALL_LIBDIR)
  set(ELFF_INSTALL_LIBDIR "${CMAKE_MPI_INSTALL_LIBDIR}")
else()
  set(ELFF_INSTALL_LIBDIR "${CMAKE_INSTALL_LIBDIR}")
endif()

if(DEFINED CMAKE_MPI_INSTALL_BINDIR AND CMAKE_MPI_INSTALL_BINDIR)
  set(ELFF_INSTALL_BINDIR "${CMAKE_MPI_INSTALL_BINDIR}")
else()
  set(ELFF_INSTALL_BINDIR "${CMAKE_INSTALL_BINDIR}")
endif()

if(DEFINED CMAKE_MPI_INSTALL_INCLUDEDIR AND CMAKE_MPI_INSTALL_INCLUDEDIR)
  set(ELFF_INSTALL_INCLUDEDIR "${CMAKE_MPI_INSTALL_INCLUDEDIR}")
else()
  set(ELFF_INSTALL_INCLUDEDIR "${CMAKE_INSTALL_INCLUDEDIR}")
endif()

#
# Build and install third-party dependencies as part of CMake build
#

add_subdirectory(third_party)  

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
set(ELFF_SOURCES "")
set(ELFF_BUILD_INCLUDE_DIR "${CMAKE_CURRENT_BINARY_DIR}/include")
set(ELFF_PUBLIC_INCLUDE_ROOT "${ELFF_BUILD_INCLUDE_DIR}/elff")

set(ELFF_DIRS 
  c/models/beam
  c/models/ibm
  config 
  general 
  io/C
  io/CXX 
  models/beam
  models/ibm 
  mpi 
  fem
)

foreach(DIR IN LISTS ELFF_DIRS)
  file(GLOB ELFF_DIR_SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/${DIR}/*.cpp")
  list(SORT ELFF_DIR_SOURCES)
  list(APPEND ELFF_SOURCES ${ELFF_DIR_SOURCES})
endforeach()

file(GLOB_RECURSE ELFF_PUBLIC_HEADERS
  CONFIGURE_DEPENDS
  RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/c/*.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/c/*.hpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/config/*.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/config/*.hpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/fem/*.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/fem/*.hpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/general/*.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/general/*.hpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/io/*.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/io/*.hpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/models/*.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/models/*.hpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/mpi/*.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/mpi/*.hpp"
)
list(SORT ELFF_PUBLIC_HEADERS)

foreach(HEADER IN LISTS ELFF_PUBLIC_HEADERS)
  get_filename_component(ELFF_HEADER_DIR
    "${ELFF_PUBLIC_INCLUDE_ROOT}/${HEADER}"
    DIRECTORY
  )
  file(MAKE_DIRECTORY "${ELFF_HEADER_DIR}")
  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/${HEADER}"
    "${ELFF_PUBLIC_INCLUDE_ROOT}/${HEADER}"
    COPYONLY
  )
endforeach()

# 
# Configure config/_config.hpp which is included in config/config.hpp
#
configure_file(
  "${PROJECT_SOURCE_DIR}/config/cmake/config.hpp.in"
  "${ELFF_PUBLIC_INCLUDE_ROOT}/config/_config.hpp"
)

#
# Declare shared (.so) library
#
add_library(ELFF SHARED ${ELFF_SOURCES})
add_library(ELFF::ELFF ALIAS ELFF) 

target_include_directories(ELFF
  PUBLIC
    $<BUILD_INTERFACE:${ELFF_BUILD_INCLUDE_DIR}>
    $<INSTALL_INTERFACE:${ELFF_INSTALL_INCLUDEDIR}>
)

set_target_properties(ELFF PROPERTIES
  OUTPUT_NAME "ELFF"
)

set_target_properties(ELFF PROPERTIES
  VERSION ${PROJECT_VERSION}
  SOVERSION ${PROJECT_VERSION_MAJOR}
)

#
# OpenMP
#
if (ELFF_USE_OPENMP)
  find_package(OpenMP COMPONENTS CXX REQUIRED)
  target_link_libraries(ELFF PUBLIC OpenMP::OpenMP_CXX)
  target_compile_definitions(ELFF PUBLIC ELFF_WITH_OPENMP=1)
endif()

#
# MPI
#
if (ELFF_USE_MPI)
  find_package(MPI COMPONENTS C CXX REQUIRED)
  target_link_libraries(ELFF PUBLIC MPI::MPI_CXX)
endif()

#
# Eigen3
#
if (ELFF_USE_EIGEN3)
  if (ELFF_BUILD_THIRD_PARTY_EIGEN3)
    # add_dependencies(ELFF eigen_ep)
    target_link_libraries(ELFF PUBLIC Eigen3::Eigen)
  else()
    find_package(Eigen3 REQUIRED)
    target_link_libraries(ELFF PUBLIC Eigen3::Eigen)
  endif()
endif()

#
# HDF5
#
if (ELFF_USE_HDF5)
  if (ELFF_USE_MPI)
    set(HDF5_PREFER_PARALLEL TRUE)
  endif()
  if (ELFF_BUILD_THIRD_PARTY_HDF5)
  else()
    find_package(HDF5 COMPONENTS C HL REQUIRED)
    target_link_libraries(ELFF PUBLIC
      hdf5::hdf5
      hdf5::hdf5_hl
    )
  endif()
endif()


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
          ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile
          WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
          COMMENT "Generating API documentation with Doxygen"
          VERBATIM)
  endif()
endif()

#
# Installation
#


# ------------------------------------------------------------------------------
# Handle optional MPI install prefix: This is entirely used for Fedora/RHEL 
# packaging, where MPI libraries are installed in a non-standard location
# ------------------------------------------------------------------------------

include(CMakePackageConfigHelpers)

# ---- 1. Install the shared library ------------------------------------------------
install(TARGETS ELFF 
  EXPORT ELFFTargets
  LIBRARY DESTINATION ${ELFF_INSTALL_LIBDIR}      
  COMPONENT ELFF_Runtime
  RUNTIME DESTINATION ${ELFF_INSTALL_BINDIR}      # libelff.so
  COMPONENT ELFF_Runtime
  INCLUDES DESTINATION ${ELFF_INSTALL_INCLUDEDIR} # for consumers
)

# ---- 2. Install all public headers -------------------------------------------------
install(DIRECTORY "${ELFF_BUILD_INCLUDE_DIR}/"
  DESTINATION ${ELFF_INSTALL_INCLUDEDIR}
  COMPONENT ELFF_Development
  FILES_MATCHING
    PATTERN "*.h"
    PATTERN "*.hpp"
)


# ---- 4. Generate CMake package configuration --------------------------------------
# This makes the project discoverable via find_package(ELFF CONFIG)
set(ELFF_INSTALL_CMAKEDIR ${ELFF_INSTALL_LIBDIR}/cmake/ELFF)

configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/config/cmake/ELFFConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/ELFFConfig.cmake
  INSTALL_DESTINATION "${ELFF_INSTALL_CMAKEDIR}"
)

write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/ELFFConfigVersion.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
)


# ---- 5. Install the package + targets ---------------------------------------------
install(EXPORT ELFFTargets
  FILE ELFFTargets.cmake
  NAMESPACE ELFF::
  DESTINATION ${ELFF_INSTALL_CMAKEDIR}
  COMPONENT ELFF_Devel
)

install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/ELFFConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/ELFFConfigVersion.cmake
  DESTINATION ${ELFF_INSTALL_CMAKEDIR}
  COMPONENT ELFF_Devel
)

# ---- 6. Export from the build tree for development -------------------------------
export(EXPORT ELFFTargets
  FILE ${CMAKE_CURRENT_BINARY_DIR}/ELFFTargets.cmake
  NAMESPACE ELFF::
)
