
include(ExternalProject)
include(GNUInstallDirs)
 
if (ELFF_USE_SUBMODULES)
  ExternalProject_Add(basilisk_ep
    UPDATE_COMMAND ""
    SOURCE_DIR "${CMAKE_SOURCE_DIR}/third_party/basilisk/basilisk_submodule"
    BINARY_DIR  "${CMAKE_BINARY_DIR}/_deps/basilisk-build"
    INSTALL_DIR "${ELFF_BASILISK_PREFIX}"
    CMAKE_ARGS
        -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
        -DCMAKE_INSTALL_BINDIR=${ELFF_INSTALL_BINDIR}
        -DCMAKE_INSTALL_LIBDIR=${ELFF_INSTALL_LIBDIR}
        -DCMAKE_INSTALL_INCLUDEDIR=${ELFF_INSTALL_INCLUDEDIR}
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
        -DCMAKE_INSTALL_BINDIR=${ELFF_INSTALL_BINDIR}
        -DCMAKE_INSTALL_LIBDIR=${ELFF_INSTALL_LIBDIR}
        -DCMAKE_INSTALL_INCLUDEDIR=${ELFF_INSTALL_INCLUDEDIR}
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
