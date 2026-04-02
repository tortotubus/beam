
include(ExternalProject)
include(GNUInstallDirs)

if (ELFF_USE_EB)
  set(ELFF_AMREX_EB ON)
else()
  set(ELFF_AMREX_EB OFF)
endif()

if (ELFF_USE_HDF5)
  set(ELFF_AMREX_HDF5 ON)
else()
  set(ELFF_AMREX_HDF5 OFF)
endif()

set(ELFF_AMREX_HDF5_ARGS
  -DAMReX_HDF5=${ELFF_AMREX_HDF5}
)
if (ELFF_USE_HDF5 AND ELFF_BUILD_THIRD_PARTY_HDF5)
  list(APPEND ELFF_AMREX_HDF5_ARGS
    -DHDF5_ROOT=${ELFF_HDF5_PREFIX}
  )
endif()
 
if (ELFF_USE_SUBMODULES)
  ExternalProject_Add(amrex_ep
    UPDATE_COMMAND ""
    SOURCE_DIR "${CMAKE_SOURCE_DIR}/third_party/amrex/amrex_submodule"
    BINARY_DIR  "${CMAKE_BINARY_DIR}/_deps/amrex-build"
    INSTALL_DIR "${ELFF_AMREX_PREFIX}"
    CMAKE_ARGS
        -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
        -DAMReX_EB=${ELFF_AMREX_EB}
        -DAMReX_PROBINIT=OFF
        ${ELFF_AMREX_HDF5_ARGS}
    UPDATE_DISCONNECTED 1
    BUILD_ALWAYS OFF
  )
else()
  ExternalProject_Add(amrex_ep
    GIT_REPOSITORY https://github.com/AMReX-Codes/amrex.git
    GIT_TAG development
    
    SOURCE_DIR  "${CMAKE_BINARY_DIR}/_deps/amrex-src"
    BINARY_DIR  "${CMAKE_BINARY_DIR}/_deps/amrex-build"
    INSTALL_DIR "${ELFF_AMREX_PREFIX}"
    CMAKE_ARGS
        -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
        -DAMReX_EB=${ELFF_AMREX_EB}
        -DAMReX_PROBINIT=OFF
        ${ELFF_AMREX_HDF5_ARGS}
    UPDATE_DISCONNECTED 1
    BUILD_ALWAYS OFF
  )
endif()

file(MAKE_DIRECTORY "${ELFF_AMREX_PREFIX}/bin")


list(PREPEND CMAKE_PREFIX_PATH "${ELFF_AMREX_PREFIX}")
