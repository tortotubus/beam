cmake_minimum_required(VERSION 3.25)

project(ELFFFluid LANGUAGES C CXX)

# Keep the solver at C++17; newer AMReX versions propagate their own higher
# language requirement through AMReX::amrex.
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(ELFF_FLUID_BUNDLED_AMREX "Build the checked-out AMReX and Hydro submodules" ON)

if(ELFF_FLUID_BUNDLED_AMREX)
  set(AMReX_SPACEDIM 3 CACHE STRING "AMReX spatial dimension")
  option(AMReX_MPI "Enable MPI" OFF)
  option(AMReX_OMP "Enable OpenMP" OFF)
  option(AMReX_EB "Enable embedded boundaries" OFF)
  option(AMReX_PARTICLES "Enable AMReX particles" OFF)
  set(AMReX_PRECISION DOUBLE CACHE STRING "AMReX floating-point precision")
  add_subdirectory(../third_party/amrex/amrex_submodule amrex)

  set(HYDRO_SPACEDIM "${AMReX_SPACEDIM}" CACHE STRING "Hydro spatial dimension")
  set(HYDRO_PRECISION "${AMReX_PRECISION}" CACHE STRING "Hydro precision")
  set(HYDRO_MPI "${AMReX_MPI}" CACHE BOOL "Hydro MPI support")
  set(HYDRO_OMP "${AMReX_OMP}" CACHE BOOL "Hydro OpenMP support")
  set(HYDRO_GPU_BACKEND "${AMReX_GPU_BACKEND}" CACHE STRING "Hydro GPU backend")
  option(HYDRO_EB "Enable Hydro embedded boundaries" OFF)
  option(HYDRO_PROJECTIONS "Enable Hydro projections" ON)
  add_subdirectory(../third_party/amrex_hydro/amrex_hydro_submodule amrex-hydro)
else()
  find_package(AMReX CONFIG REQUIRED)
  find_package(AMReX-Hydro CONFIG REQUIRED)
endif()

file(GLOB FLUID_SOURCES CONFIGURE_DEPENDS src/*.cpp src/dlm/*.cpp)
add_library(elff_fluid STATIC ${FLUID_SOURCES})
add_library(ELFF::fluid ALIAS elff_fluid)
set_target_properties(elff_fluid PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(elff_fluid PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_compile_features(elff_fluid PUBLIC cxx_std_17)
target_link_libraries(elff_fluid PUBLIC AMReX::amrex AMReX-Hydro::amrex_hydro_api)

add_executable(elff_fluid_solver app/main.cpp)
target_link_libraries(elff_fluid_solver PRIVATE ELFF::fluid)
