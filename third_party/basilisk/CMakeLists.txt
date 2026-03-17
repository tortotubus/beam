
include(ExternalProject)
include(GNUInstallDirs)

message(STATUS "ELFF_BASILISK_PREFIX=${ELFF_BASILISK_PREFIX}")


if (ELFF_USE_SUBMODULES)
  ExternalProject_Add(basilisk_ep
    UPDATE_COMMAND ""
    SOURCE_DIR "${CMAKE_SOURCE_DIR}/third_party/basilisk/basilisk_submodule"
    BINARY_DIR  "${CMAKE_BINARY_DIR}/_deps/basilisk-build"
    INSTALL_DIR "${ELFF_BASILISK_PREFIX}"
    CMAKE_ARGS
        -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>      
    UPDATE_DISCONNECTED 1
    BUILD_ALWAYS OFF
  )
else()
  ExternalProject_Add(basilisk_ep
    GIT_REPOSITORY https://github.com/tortotubus/pacific-basilisk.git

    SOURCE_DIR  "${CMAKE_BINARY_DIR}/_deps/basilisk-src"
    BINARY_DIR  "${CMAKE_BINARY_DIR}/_deps/basilisk-build"
    INSTALL_DIR "${ELFF_BASILISK_PREFIX}"
    CMAKE_ARGS
        -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>        
    UPDATE_DISCONNECTED 1
    BUILD_ALWAYS OFF
  )
endif()

file(MAKE_DIRECTORY "${ELFF_BASILISK_PREFIX}/bin")

add_executable(basilisk::qcc IMPORTED GLOBAL)
set_target_properties(basilisk::qcc PROPERTIES
  IMPORTED_LOCATION "${ELFF_BASILISK_PREFIX}/bin/qcc"
)
add_dependencies(basilisk::qcc basilisk_ep)


list(PREPEND CMAKE_PREFIX_PATH "${ELFF_BASILISK_PREFIX}")
