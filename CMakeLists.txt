cmake_minimum_required(VERSION 3.10...3.27)

# Enable both C and C++
project(beam VERSION 1.0 LANGUAGES C CXX)

# Set CMake policies
cmake_policy(SET CMP0167 NEW)  # Use FindBoost from CMake 3.27

# —————————————————————————————————————————————
# Options & flags
# —————————————————————————————————————————————
# Build shared by default; set -DBUILD_SHARED_LIBS=OFF to get a static .a only
option(BUILD_SHARED_LIBS "Build shared libraries" ON)

# Set C++17 for if constexpr support
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# PIC for static libs so they can be linked into shared objects
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Include path
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

# Find Eigen Package
find_package(Eigen3 3.3 REQUIRED)

# Find VTK package with additional components
find_package(VTK 9.0 COMPONENTS 
    CommonCore
    CommonDataModel
    IOXMLParser
    IOXML
    RenderingCore
    RenderingOpenGL2
    InteractionStyle
    RenderingContextOpenGL2
    InteractionWidgets
    REQUIRED
)

# Find MFEM package
find_package(MFEM REQUIRED)

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
    ${VTK_LIBRARIES}
    ${MFEM_LIBRARIES}
    Eigen3::Eigen
)

target_include_directories(beam PUBLIC ${VTK_INCLUDE_DIRS} ${MFEM_INCLUDE_DIRS})

# VTK requires this when building with static libraries
if(VTK_ENABLE_WRAPPING)
    vtk_module_autoinit(
        TARGETS beam
        MODULES ${VTK_LIBRARIES}
    )
endif()

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

# Add executable
add_executable(beam_demo src/main.cpp)
target_link_libraries(beam_demo 
    PRIVATE 
    beam
    ${VTK_LIBRARIES}
    ${MFEM_LIBRARIES}
)

# VTK module initialization for the demo executable
if(VTK_ENABLE_WRAPPING)
    vtk_module_autoinit(
        TARGETS beam_demo
        MODULES ${VTK_LIBRARIES}
    )
endif()

# Add Doxygen documentation
find_package(Doxygen)
if(DOXYGEN_FOUND)
    add_custom_target(doc
        ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        COMMENT "Generating API documentation with Doxygen"
        VERBATIM)
endif()
