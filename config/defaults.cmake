
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
      "Build type: Debug, Release, RelWithDebInfo, or MinSizeRel." FORCE)
endif()

option(BUILD_SHARED_LIBS "Enable shared library build of BEAM" ON)
set(ELFF_PRECISION "double" CACHE STRING "Floating-point precision to use: single, or double")
option(ELFF_USE_EXCEPTIONS "Enable the use of exceptions" OFF)
option(ELFF_USE_MPI "Enable the use of MPI" ON)
option(ELFF_USE_OPENMP "Enable the use of OpenMP" OFF)
option(ELFF_USE_HDF5 "Enable the use of HDF5" ON)
option(ELFF_USE_EIGEN3 "Enable the use of Eigen3" ON)
option(ELFF_USE_BASILISK "Enable the basilisk subproject" ON)

option(ELFF_BUILD_DOC "Build documentation" OFF)
option(ELFF_BUILD_TESTS "Build unit & integration tests" OFF)

option(ELFF_USE_SUBMODULES "Use git submodules for third-party dependencies where available" OFF)
option(ELFF_BUILD_THIRD_PARTY_HDF5 "Build HDF5" OFF)
option(ELFF_BUILD_THIRD_PARTY_BASILISK "Build Basilisk" ON)
option(ELFF_BUILD_THIRD_PARTY_EIGEN3 "Build Eigen" OFF)