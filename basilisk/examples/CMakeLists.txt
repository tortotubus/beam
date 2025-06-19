
set(CMAKE_VERBOSE_MAKEFILE ON)

# Enable C language
enable_language(C)

# Locate QCC on the user’s PATH (or via -DQCC_EXECUTABLE)
find_program(
  QCC_EXECUTABLE
  NAMES qcc
  HINTS ENV PATH
)
if (NOT QCC_EXECUTABLE)
  message(FATAL_ERROR
    "Could not find 'qcc'.\n"
    "Please install QCC or set PATH so that 'qcc' is on your shell's path,\n"
    "or pass -DQCC_EXECUTABLE=/full/path/to/qcc to CMake.")
endif()

message(STATUS "Using QCC compiler at ${QCC_EXECUTABLE}")

# -------------------------------------------------------------------
# 1) Copy interface headers at build time
# -------------------------------------------------------------------
# We'll make a custom target that copies the entire directory
file(GLOB BASILISK_INTERFACE_HDRS
     "${CMAKE_SOURCE_DIR}/basilisk/interface/*.h"
)

file(MAKE_DIRECTORY
     "${CMAKE_BINARY_DIR}/basilisk/interface"
)

file(COPY ${BASILISK_INTERFACE_HDRS}
     DESTINATION "${CMAKE_BINARY_DIR}/basilisk/interface"
)

add_custom_target(copy_interface_headers ALL
  COMMAND ${CMAKE_COMMAND} -E make_directory
          "${CMAKE_BINARY_DIR}/basilisk/interface"
  # only copy the .h files we globbed above:
  COMMAND ${CMAKE_COMMAND} -E copy 
          ${BASILISK_INTERFACE_HDRS} 
          "${CMAKE_BINARY_DIR}/basilisk/interface/"
  COMMENT "Copying basilisk/interface/*.h → build/basilisk/interface/"
  VERBATIM
)

# -------------------------------------------------------------------
# 2) Copy .c examples at build time
# -------------------------------------------------------------------
file(GLOB BASILISK_C_EXAMPLE_SRCS
     "${CMAKE_SOURCE_DIR}/basilisk/examples/*.c"
)
file(MAKE_DIRECTORY
     "${CMAKE_BINARY_DIR}/basilisk/examples"
)
file(COPY ${BASILISK_C_EXAMPLE_SRCS}
     DESTINATION "${CMAKE_BINARY_DIR}/basilisk/examples"
)
add_custom_target(copy_examples ALL
  COMMAND ${CMAKE_COMMAND} -E make_directory
          "${CMAKE_BINARY_DIR}/basilisk/examples"
  # only copy the .c files we globbed above:
  COMMAND ${CMAKE_COMMAND} -E copy
          ${BASILISK_C_EXAMPLE_SRCS}
          "${CMAKE_BINARY_DIR}/basilisk/examples/"
  COMMENT "Copying basilisk/examples/*.c → build/basilisk/examples/"
  VERBATIM
)

# -------------------------------------------------------------------
# 3) Now glob the copied .c files and build each one
# -------------------------------------------------------------------
file(GLOB BASILISK_C_EXAMPLE_SRCS_CPY
     "${CMAKE_BINARY_DIR}/basilisk/examples/*.c"
)

file(GLOB BASILISK_C_EXAMPLE_SRCS_ORIG
     "${CMAKE_SOURCE_DIR}/basilisk/examples/*.c"
)

foreach(example_basilisk_source ${BASILISK_C_EXAMPLE_SRCS_ORIG})
  get_filename_component(example_name ${example_basilisk_source} NAME_WE)

  # Tell CMake how to generate that file
  add_custom_command(
    OUTPUT ${example_name}
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/basilisk/examples/"
    COMMAND ${QCC_EXECUTABLE} "${example_name}.c" -o "${example_name}" -L${CMAKE_BINARY_DIR} -lm -lbeam -Wl,-rpath,'${CMAKE_BINARY_DIR}'
    DEPENDS "${example_basilisk_source}"
  )

  add_custom_target(${example_name}_exe ALL DEPENDS ${example_name} beam copy_interface_headers copy_examples)
  
endforeach()
