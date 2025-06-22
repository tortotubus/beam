
set(CMAKE_VERBOSE_MAKEFILE ON)

# Enable C language
enable_language(C)

set(MY_MAKEFILE "${CMAKE_BINARY_DIR}/basilisk/examples/Makefile")
file(WRITE "${MY_MAKEFILE}"
  "include \$(BASILISK)/Makefile.defs\n"
  "LIBS += -I${CMAKE_SOURCE_DIR}/basilisk/ -lhdf5 -lhdf5_hl "
    "-L${CMAKE_BINARY_DIR} -Wl,-rpath,'${CMAKE_BINARY_DIR}' -lbeam\n"
)
file(GLOB BASILISK_C_EXAMPLE_SRCS "${CMAKE_SOURCE_DIR}/basilisk/examples/*.c")
file(GLOB_RECURSE BASILISK_LIBRARY_HDRS "${CMAKE_SOURCE_DIR}/basilisk/library/*.h")

foreach(example_basilisk_source ${BASILISK_C_EXAMPLE_SRCS})
  get_filename_component(example_name ${example_basilisk_source} NAME_WE)

  # Tell CMake how to generate that file
  add_custom_command(
    OUTPUT "${CMAKE_BINARY_DIR}/basilisk/examples/${example_name}"
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/basilisk/examples"
    # 1) ensure the build dir exists
    COMMAND ${CMAKE_COMMAND} -E make_directory
            "${CMAKE_BINARY_DIR}/basilisk/examples"
    # 2) copy the .c file over
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${CMAKE_SOURCE_DIR}/basilisk/examples/${example_name}.c"
            "${CMAKE_BINARY_DIR}/basilisk/examples/${example_name}.c"
    # 3) compile
    COMMAND ${QCC_EXECUTABLE}
      -autolink
      -nolineno
      -s
      #-dimensions=dims
      -non-finite
      -Wall -D_FORTIFY_SOURCE=2
      "${example_name}.c" 
      -o "${example_name}"
      -I"${CMAKE_SOURCE_DIR}/basilisk/"
      -g -pipe 
      -lhdf5 -lhdf5_hl
      -L"${CMAKE_BINARY_DIR}"
      -Wl,-rpath,'${CMAKE_BINARY_DIR}'
      -lbeam 
      -lm 
    # 4) delete the .c file
    # COMMAND ${CMAKE_COMMAND} -E rm  "${CMAKE_BINARY_DIR}/basilisk/examples/${example_name}.c"
    DEPENDS
      ${BASILISK_LIBRARY_HDRS}
      ${example_basilisk_source}
      beam
  )

  add_custom_target(${example_name}.build ALL
    DEPENDS "${CMAKE_BINARY_DIR}/basilisk/examples/${example_name}"
  )

endforeach()
