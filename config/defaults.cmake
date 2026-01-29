
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
      "Build type: Debug, Release, RelWithDebInfo, or MinSizeRel." FORCE)
endif()

option(BUILD_SHARED_LIBS "Enable shared library build of BEAM" OFF)

set(ELFF_PRECISION "double" CACHE STRING "Floating-point precision to use: single, or double")

option(ELFF_USE_EXCEPTIONS "Enable the use of exceptions" OFF)

