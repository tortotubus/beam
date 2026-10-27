include(ExternalProject)
include(GNUInstallDirs)

if (ELFF_USE_SUBMODULES)
  ExternalProject_Add(amrex_ep
    UPDATE_COMMAND ""
    SOURCE_DIR "${CMAKE_SOURCE_DIR}/third_party/amrex/amrex_submodule" 
    BINARY_DIR  "${CMAKE_BINARY_DIR}/_deps/amrex-build"
    INSTALL_DIR "${ELFF_AMREX_PREFIX}"
    CMAKE_ARGS
      -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
      -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        
    UPDATE_DISCONNECTED 1
    BUILD_ALWAYS OFF
  )
else()
  ExternalProject_Add(amrex_ep
    GIT_REPOSITORY https://github.com/AMReX-Codes/amrex.git
    GIT_TAG        development

    SOURCE_DIR  "${CMAKE_BINARY_DIR}/_deps/amrex-src"
    BINARY_DIR  "${CMAKE_BINARY_DIR}/_deps/amrex-build"
    INSTALL_DIR "${ELFF_AMREX_PREFIX}"
  
    CMAKE_ARGS
      -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
      -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
      
    UPDATE_DISCONNECTED 1
    BUILD_ALWAYS OFF
  )
endif()

file(MAKE_DIRECTORY "${ELFF_AMREX_PREFIX}/include")
file(MAKE_DIRECTORY "${ELFF_AMREX_PREFIX}/${CMAKE_INSTALL_LIBDIR}")


list(PREPEND CMAKE_PREFIX_PATH "${ELFF_AMREX_PREFIX}")
