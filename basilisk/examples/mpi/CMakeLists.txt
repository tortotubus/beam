
set(CMAKE_VERBOSE_MAKEFILE ON)

# Enable C language
enable_language(C)


find_package(MPI REQUIRED)
find_package(HDF5 REQUIRED COMPONENTS C HL)

if(NOT HDF5_IS_PARALLEL)
  message(FATAL_ERROR "Need parallel HDF5 (built with MPI-IO).")
endif()

file(GLOB BASILISK_C_EXAMPLE_SRCS "${CMAKE_SOURCE_DIR}/basilisk/examples/mpi/*.c")
file(GLOB_RECURSE BASILISK_LIBRARY_HDRS "${CMAKE_SOURCE_DIR}/basilisk/library/*.h")

foreach(example_basilisk_source ${BASILISK_C_EXAMPLE_SRCS})
  get_filename_component(example_name ${example_basilisk_source} NAME_WE)

  add_custom_command(
    OUTPUT "${CMAKE_BINARY_DIR}/basilisk/examples/mpi/_${example_name}.c"
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/basilisk/examples/mpi"
    # 1) ensure the build dir exists
    COMMAND ${CMAKE_COMMAND} -E make_directory
            "${CMAKE_BINARY_DIR}/basilisk/examples"
    # 2) copy the .c file over
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${CMAKE_SOURCE_DIR}/basilisk/examples/mpi/${example_name}.c"
            "${CMAKE_BINARY_DIR}/basilisk/examples/mpi/${example_name}.c"
    # 3) compile
    COMMAND ${QCC_EXECUTABLE}
      -D_MPI=1
      "${example_name}.c" 
      -I"${CMAKE_SOURCE_DIR}/basilisk/"
      # -DTRACE=1
      -source
    DEPENDS
      basilisk #source
      ${BASILISK_LIBRARY_HDRS}
      ${example_basilisk_source}
  )


  add_executable(${example_name} "_${example_name}.c")

  set_target_properties(${example_name} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/basilisk/examples/mpi"
  )


  # 4. Include your headers
  target_include_directories(${example_name} PRIVATE
    "${CMAKE_SOURCE_DIR}/basilisk"
  )

  target_compile_options(${example_name} PRIVATE
    -Wall
    -D_FORTIFY_SOURCE=2
    -g
    -pipe
  )

  target_link_libraries(${example_name}
    PRIVATE
      beam
      MPI::MPI_C
      hdf5::hdf5
      hdf5::hdf5_hl
      m
  )
  
  set_target_properties(${example_name} PROPERTIES
    BUILD_RPATH "${CMAKE_BINARY_DIR}"
  )

  # add_custom_target(mpirun_8_${example_name}
  #   COMMAND ${MPIEXEC_EXECUTABLE}
  #           ${MPIEXEC_NUMPROC_FLAG} 8
  #           ${MPIEXEC_PREFLAGS}
  #           $<TARGET_FILE:${example_name}>
  #           ${MPIEXEC_POSTFLAGS}
  #   DEPENDS ${example_name}
  #   WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/basilisk/examples/mpi"
  #   USES_TERMINAL
  #   COMMENT "Running ${example_name} with 8 MPI ranks"
  # )

endforeach()
