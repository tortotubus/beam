option(BUILD_BASILISK_INTERFACE "Build basilisk interface" ON)

if (BUILD_BASILISK_INTERFACE)
  add_subdirectory(interface)
endif()

#
#
#

option(BUILD_BASILISK_EXAMPLES "Build basilisk example programs" ON)

if(BUILD_BASILISK_EXAMPLES)
  find_package(basilisk CONFIG REQUIRED) 
  include(elff_basilisk_add_executable)
  add_subdirectory(examples/serial)
endif()

option(BUILD_BASILISK_TESTS "Build basilisk example programs" ON)

if(BUILD_BASILISK_TESTS)
  find_package(basilisk CONFIG REQUIRED) 
  include(elff_basilisk_add_executable)
  add_subdirectory(tests)
endif()


set(SOURCES ${SOURCES} ${SRCS} PARENT_SCOPE)
set(HEADERS ${HEADERS} ${HDRS} PARENT_SCOPE)
