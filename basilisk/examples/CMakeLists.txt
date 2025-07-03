
set(CMAKE_VERBOSE_MAKEFILE ON)

# Enable C language
enable_language(C)

# Find HDF5 package
find_package(HDF5)
# 1) Find MPI
find_package(MPI REQUIRED)

file(GLOB BASILISK_C_EXAMPLE_SRCS "${CMAKE_SOURCE_DIR}/basilisk/examples/*.c")
file(GLOB_RECURSE BASILISK_LIBRARY_HDRS "${CMAKE_SOURCE_DIR}/basilisk/library/*.h")

foreach(example_basilisk_source ${BASILISK_C_EXAMPLE_SRCS})
  get_filename_component(example_name ${example_basilisk_source} NAME_WE)

  add_custom_command(
    OUTPUT "${CMAKE_BINARY_DIR}/basilisk/examples/_${example_name}.c"
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/basilisk/examples"
    # 1) ensure the build dir exists
    COMMAND ${CMAKE_COMMAND} -E make_directory
            "${CMAKE_BINARY_DIR}/basilisk/examples"
    # 2) copy the .c file over
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${CMAKE_SOURCE_DIR}/basilisk/examples/${example_name}.c"
            "${CMAKE_BINARY_DIR}/basilisk/examples/${example_name}.c"
    # 3) compile
    COMMAND $<TARGET_FILE:qcc-source>
      "${example_name}.c" 
      -I"${CMAKE_SOURCE_DIR}/basilisk/"
      -I "${MPI_INCLUDE_PATHS}"
      -I "${HDF5_INCLUDE_DIRS}"
      -source
    DEPENDS
      qcc-source
      ${BASILISK_LIBRARY_HDRS}
      ${example_basilisk_source}
  )


  add_executable(${example_name} "_${example_name}.c")

  set_target_properties(${example_name} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/basilisk/examples"
  )


  # 4. Include your headers
  target_include_directories(${example_name} PRIVATE
    "${CMAKE_SOURCE_DIR}/basilisk"
    ${MPI_INCLUDE_PATHS}
    ${HDF5_INCLUDE_DIRS}
  )

  target_compile_options(${example_name} PRIVATE
    #-autolink
    -Wall
    -D_FORTIFY_SOURCE=2
    -g
    -pipe
  )

  target_link_libraries(${example_name}
    PRIVATE
      beam
      ${HDF5_HL_LIBRARIES}
      ${HDF5_LIBRARIES}
      # hdf5::hdf5
      # hdf5::hdf5_hl
      m
  )
  
  set_target_properties(${example_name} PROPERTIES
    BUILD_RPATH "${CMAKE_BINARY_DIR}"
  )

endforeach()
