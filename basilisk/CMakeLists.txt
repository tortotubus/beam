option(BUILD_BASILISK_INTERFACE "Build basilisk interface" ON)
if (BUILD_BASILISK_INTERFACE)
  add_subdirectory(interface)
endif()

option(BUILD_BASILISK_EXAMPLES "Build basilisk example programs" ON)
if(BUILD_BASILISK_EXAMPLES)
  add_subdirectory(examples)
endif()

option(BUILD_BASILISK_QCC "Build basilisk compiler" OFF)

if(BUILD_BASILISK_QCC)
  # # find GNU make
  # find_program(MAKE_PROG NAMES make REQUIRED
  #   DOC "The GNU Make executable")

  # # real source & out-of-tree copy
  # set(BASILISK_SRC   ${CMAKE_SOURCE_DIR}/basilisk/basilisk/src)
  # set(BASILISK_BUILD ${CMAKE_BINARY_DIR}/external/basilisk/src)

  # add_custom_target(build_basilisk ALL
  #   COMMENT "Fresh-copy & build Basilisk out-of-tree"
  #   # 1) remove any previous copy
  #   COMMAND ${CMAKE_COMMAND} -E remove_directory ${BASILISK_BUILD}
  #   # 2) recreate the directory
  #   COMMAND ${CMAKE_COMMAND} -E make_directory  ${BASILISK_BUILD}
  #   # 3) copy pristine sources
  #   COMMAND ${CMAKE_COMMAND} -E copy_directory ${BASILISK_SRC} ${BASILISK_BUILD}
  #   # 4) symlink config→config.gcc
  #   COMMAND ${CMAKE_COMMAND} -E remove -f ${BASILISK_BUILD}/config
  #   COMMAND ${CMAKE_COMMAND} -E create_symlink config.gcc ${BASILISK_BUILD}/config
  #   # 5) invoke GNU make in the copy
  #   COMMAND ${MAKE_PROG} -C ${BASILISK_BUILD}
  #   VERBATIM
  # )
endif()

set(SOURCES ${SOURCES} ${SRCS} PARENT_SCOPE)
set(HEADERS ${HEADERS} ${HDRS} PARENT_SCOPE)
