cmake_minimum_required(VERSION 3.18)

if (CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  project(IAMR LANGUAGES C CXX)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

get_filename_component(IAMR_REPO_ROOT "${CMAKE_CURRENT_LIST_DIR}/.." ABSOLUTE)
set(IAMR_LOCAL_PREFIX "${IAMR_REPO_ROOT}/install/release")

if (EXISTS "${IAMR_LOCAL_PREFIX}")
  list(PREPEND CMAKE_PREFIX_PATH "${IAMR_LOCAL_PREFIX}")
endif()

if (NOT TARGET AMReX::amrex)
  find_package(AMReX CONFIG REQUIRED 3D DOUBLE MPI)
endif()

if (NOT TARGET AMReX-Hydro::amrex_hydro_api)
  find_package(AMReX-Hydro CONFIG REQUIRED)
endif()

set(IAMR_SOURCES
  src/main.cpp
  src/Diffusion.cpp
  src/FluxBoxes.cpp
  src/MacProj.cpp
  src/NavierStokes.cpp
  src/NavierStokesBase.cpp
  src/NS_average.cpp
  src/NS_derive.cpp
  src/NS_error.cpp
  src/NS_getForce.cpp
  src/NS_init_eb2.cpp
  src/NS_LES.cpp
  src/NS_setup.cpp
  src/NS_util.cpp
  src/NSBld.cpp
  src/OutFlowBC.cpp
  src/Projection.cpp
  src/SyncRegister.cpp
  src/Utilities/Utilities.cpp
  src/prob/prob_init.cpp
)

add_executable(iamr)
target_sources(iamr PRIVATE ${IAMR_SOURCES})

target_include_directories(iamr
  PRIVATE
    "${CMAKE_CURRENT_LIST_DIR}/src"
    "${CMAKE_CURRENT_LIST_DIR}/src/Utilities"
    "${CMAKE_CURRENT_LIST_DIR}/src/prob"
)

target_link_libraries(iamr
  PRIVATE
    AMReX::amrex
    AMReX-Hydro::amrex_hydro_api
)

if (TARGET amrex_ep)
  add_dependencies(iamr amrex_ep)
endif()

if (TARGET amrex-hydro_ep)
  add_dependencies(iamr amrex-hydro_ep)
endif()

include(AMReXBuildInfo)
generate_buildinfo(iamr "${IAMR_REPO_ROOT}")
