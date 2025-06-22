option(BUILD_BASILISK_INTERFACE "Build basilisk interface" ON)
if (BUILD_BASILISK_INTERFACE)
  add_subdirectory(interface)
endif()

option(BUILD_BASILISK_EXAMPLES "Build basilisk example programs" ON)
if(BUILD_BASILISK_EXAMPLES)
  add_subdirectory(examples)
endif()

include(ExternalProject)

ExternalProject_Add(qcc_ep
  SOURCE_DIR    ${CMAKE_SOURCE_DIR}/basilisk/basilisk
  BINARY_DIR    ${CMAKE_BINARY_DIR}/basilisk/basilisk/build
  INSTALL_DIR   ${CMAKE_BINARY_DIR}/basilisk/basilisk/build/install
  CMAKE_ARGS
    -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
)

ExternalProject_Get_Property(qcc_ep install_dir)

# Locate QCC on the user’s PATH (or via -DQCC_EXECUTABLE)
find_program(
  QCC_EXECUTABLE
  NAMES qcc
  HINTS HINTS ${install_dir}/bin
)

if (NOT QCC_EXECUTABLE)
  message(FATAL_ERROR
    "Could not find 'qcc'.\n"
    "Please install QCC or set PATH so that 'qcc' is on your shell's path,\n"
    "or pass -DQCC_EXECUTABLE=/full/path/to/qcc to CMake.")
endif()

message(STATUS "Using QCC compiler at ${QCC_EXECUTABLE}")


set(SOURCES ${SOURCES} ${SRCS} PARENT_SCOPE)
set(HEADERS ${HEADERS} ${HDRS} PARENT_SCOPE)
