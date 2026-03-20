

option(BUILD_BASILISK_EXAMPLES "Build basilisk example programs" ON)

if(BUILD_BASILISK_EXAMPLES)
  if(ELFF_USE_BASILISK)
  else()
  find_package(basilisk CONFIG REQUIRED) 
  endif()
  find_package(MPI REQUIRED)
  include(elff_basilisk_add_executable)
  include(elff_basilisk_mpi_add_executable)
  add_subdirectory(examples/serial)
  add_subdirectory(examples/mpi)
endif()


set(SOURCES ${SOURCES} ${SRCS} PARENT_SCOPE)
set(HEADERS ${HEADERS} ${HDRS} PARENT_SCOPE)
