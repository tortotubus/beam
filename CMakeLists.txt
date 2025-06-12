cmake_minimum_required(VERSION 3.5)
project(beam VERSION 1.0 LANGUAGES CXX)

# —————————————————————————————————————————————
# Options & flags
# —————————————————————————————————————————————
# Build shared by default; set -DBUILD_SHARED_LIBS=OFF to get a static .a only
option(BUILD_SHARED_LIBS "Build shared libraries" ON)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# PIC for static libs so they can be linked into shared objects
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Include path
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

# —————————————————————————————————————————————
# Locate LAPACK/BLAS/LAPACKE
# —————————————————————————————————————————————
#find_package(LAPACK REQUIRED)   # provides LAPACK_LIBRARIES
#find_package(BLAS   REQUIRED)   # provides BLAS_LIBRARIES
#find_library(LAPACKE_LIBRARY
#  NAMES lapacke
#  REQUIRED
#)

find_package(Eigen3 3.3 REQUIRED)

# —————————————————————————————————————————————
# Library target
# —————————————————————————————————————————————
set(SRCS
  src/Euler/InextensibleEulerBeam.cpp
  src/beam_c.cpp
  src/LinAlg/Matrix.cpp
  src/LinAlg/Vector.cpp
  src/MyClass.cpp
)

add_library(beam ${SRCS})

# Link in numeric libs
target_link_libraries(beam
  PUBLIC
    #${LAPACKE_LIBRARY}
    #${LAPACK_LIBRARIES}
    #${BLAS_LIBRARIES}
    Eigen3::Eigen
    autodiff:autodiff
)

# Versioning
set_target_properties(beam PROPERTIES
  VERSION     ${PROJECT_VERSION}
  SOVERSION   ${PROJECT_VERSION_MAJOR}
)

# —————————————————————————————————————————————
# Install rules
# —————————————————————————————————————————————
install(TARGETS beam
        EXPORT beam-config
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        RUNTIME DESTINATION bin)

install(DIRECTORY include/ DESTINATION include)

install(EXPORT beam-config
        NAMESPACE beam::
        DESTINATION lib/cmake/beam)


find_package(GTest REQUIRED)
enable_testing()

option(BUILD_TESTS "Build unit & integration tests" ON)
if(BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()

# —————————————————————————————————————————————
# pkg-config support (optional)
# —————————————————————————————————————————————
# find_package(PkgConfig QUIET)
# if(PKG_CONFIG_FOUND)
#   pkg_check_modules(PC_BEAM QUIET beam)
#   if(NOT PC_BEAM_FOUND)
#     # generate a simple pc file
#     configure_file(${CMAKE_CURRENT_SOURCE_DIR}/beam.pc.in
#                    ${CMAKE_CURRENT_BINARY_DIR}/beam.pc @ONLY)
#     install(FILES ${CMAKE_CURRENT_BINARY_DIR}/beam.pc
#             DESTINATION lib/pkgconfig)
#   endif()
# endif()

# Add executable
add_executable(beam_demo src/main.cpp)
target_link_libraries(beam_demo PRIVATE beam)
