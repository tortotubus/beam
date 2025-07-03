option(BUILD_BASILISK_INTERFACE "Build basilisk interface" ON)
if (BUILD_BASILISK_INTERFACE)
  add_subdirectory(interface)
endif()


# 2) (Optional) Define where you want the install tree of the external project
# set(BEAM_BASILISK_INSTALL_DIR
#     ${CMAKE_BINARY_DIR}/basilisk/basilisk/build/install
# )

# 3) Tell FetchContent how to grab the repo
#    You can also pass CMake arguments by setting <name>_CMAKE_ARGS
# set(beam_basilisk_CMAKE_ARGS
#   -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
#   -DCMAKE_INSTALL_PREFIX=${BEAM_BASILISK_INSTALL_DIR}
# )


# include(ExternalProject)

# ExternalProject_Add(qcc_ep
#   GIT_REPOSITORY git@gitlab.math.ubc.ca:colive/beam-basilisk
#   UPDATE_DISCONNECTED 1
#   SOURCE_DIR    ${CMAKE_SOURCE_DIR}/basilisk/basilisk
#   BINARY_DIR    ${CMAKE_BINARY_DIR}/basilisk/basilisk/build
#   INSTALL_DIR   ${CMAKE_BINARY_DIR}/basilisk/basilisk/build/install
#   CMAKE_ARGS
#     -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
#     -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
# )

# ExternalProject_Add(qcc_ep
#   GIT_REPOSITORY        git@gitlab.math.ubc.ca:colive/beam-basilisk
#   UPDATE_DISCONNECTED   1
#   SOURCE_DIR            ${CMAKE_SOURCE_DIR}/basilisk/basilisk
#   BINARY_DIR            ${CMAKE_BINARY_DIR}/basilisk/basilisk/build
#   INSTALL_DIR           ${CMAKE_BINARY_DIR}/basilisk/basilisk/build/install

#   CMAKE_ARGS
#     -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
#     -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
 
#   BUILD_COMMAND
#     ${CMAKE_COMMAND}
#       --build <BINARY_DIR>
#       --config ${CMAKE_BUILD_TYPE}
#       --target qcc
# )

# ExternalProject_Get_Property(qcc_ep install_dir)

# 1) Enable FetchContent
include(FetchContent)

FetchContent_Declare(
  basilisk
  GIT_REPOSITORY git@gitlab.math.ubc.ca:colive/beam-basilisk.git
)

# 4) Bring it into your build
FetchContent_MakeAvailable(basilisk)

# set(QCC_EXECUTABLE "${install_dir}/bin/qcc" CACHE FILEPATH "Path to built qcc")
# add_custom_target(ensure-qcc DEPENDS qcc_ep)

option(BUILD_BASILISK_EXAMPLES "Build basilisk example programs" ON)
if(BUILD_BASILISK_EXAMPLES)
  add_subdirectory(examples)
endif()


set(SOURCES ${SOURCES} ${SRCS} PARENT_SCOPE)
set(HEADERS ${HEADERS} ${HDRS} PARENT_SCOPE)
