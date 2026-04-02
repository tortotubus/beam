
include(ExternalProject)
include(GNUInstallDirs)

if (ELFF_USE_EB)
  set(ELFF_HYDRO_EB ON)
  set(ELFF_HYDRO_NO_EB OFF)
else()
  set(ELFF_HYDRO_EB OFF)
  set(ELFF_HYDRO_NO_EB ON)
endif()

set(ELFF_HYDRO_HDF5_ARGS "")
if (ELFF_USE_HDF5 AND ELFF_BUILD_THIRD_PARTY_HDF5)
  list(APPEND ELFF_HYDRO_HDF5_ARGS
    -DHDF5_ROOT=${ELFF_HDF5_PREFIX}
  )
endif()
 
if (ELFF_USE_SUBMODULES)
  ExternalProject_Add(amrex-hydro_ep
    UPDATE_COMMAND ""
    SOURCE_DIR "${CMAKE_SOURCE_DIR}/third_party/amrex-hydro/amrex-hydro_submodule"
    BINARY_DIR  "${CMAKE_BINARY_DIR}/_deps/amrex-hydro-build"
    INSTALL_DIR "${ELFF_AMREX_HYDRO_PREFIX}"
    CMAKE_ARGS
        -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
        -DHYDRO_EB=${ELFF_HYDRO_EB}
        -DHYDRO_NO_EB=${ELFF_HYDRO_NO_EB}
        -DAMReX_DIR=${ELFF_AMREX_PREFIX}/lib/cmake/AMReX
        ${ELFF_HYDRO_HDF5_ARGS}
    UPDATE_DISCONNECTED 1
    BUILD_ALWAYS OFF
    DEPENDS amrex_ep
  )
else()
  ExternalProject_Add(amrex-hydro_ep
    GIT_REPOSITORY https://github.com/AMReX-Codes/AMReX-Hydro.git
    GIT_TAG development

    SOURCE_DIR  "${CMAKE_BINARY_DIR}/_deps/amrex-hydro-src"
    BINARY_DIR  "${CMAKE_BINARY_DIR}/_deps/amrex-hydro-build"
    INSTALL_DIR "${ELFF_AMREX_HYDRO_PREFIX}"
    CMAKE_ARGS
        -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
        -DHYDRO_EB=${ELFF_HYDRO_EB}
        -DHYDRO_NO_EB=${ELFF_HYDRO_NO_EB}
        -DAMReX_DIR=${ELFF_AMREX_PREFIX}/lib/cmake/AMReX
        ${ELFF_HYDRO_HDF5_ARGS}
    UPDATE_DISCONNECTED 1
    BUILD_ALWAYS OFF
    DEPENDS amrex_ep
  )
endif()




list(PREPEND CMAKE_PREFIX_PATH "${ELFF_AMREX_HYDRO_PREFIX}")
