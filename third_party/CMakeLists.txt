 
if (ELFF_USE_HDF5)
    if (ELFF_BUILD_THIRD_PARTY_HDF5)
        set(ELFF_HDF5_PREFIX "${CMAKE_INSTALL_PREFIX}" CACHE PATH "HDF5 install prefix")
        add_subdirectory(hdf5)
    endif()
endif()


if (ELFF_USE_BASILISK)
    if (ELFF_BUILD_THIRD_PARTY_BASILISK)
        set(ELFF_BASILISK_PREFIX "${CMAKE_INSTALL_PREFIX}" CACHE PATH "Basilisk install prefix")
        add_subdirectory(basilisk)
    endif()
endif()


if (ELFF_USE_EIGEN3)
    if (ELFF_BUILD_THIRD_PARTY_EIGEN3)
        set(ELFF_EIGEN3_PREFIX "${CMAKE_INSTALL_PREFIX}" CACHE PATH "Eigen3 install prefix")
        add_subdirectory(eigen)
    endif()
endif()

if (ELFF_USE_AMREX)
    if (ELFF_BUILD_THIRD_PARTY_AMREX)
        set(ELFF_AMREX_PREFIX "${CMAKE_INSTALL_PREFIX}" CACHE PATH "AMReX install prefix")
        add_subdirectory(amrex)
    endif()
endif()

if (ELFF_USE_AMREX_HYDRO)
    if (ELFF_BUILD_THIRD_PARTY_AMREX_HYDRO)
        set(ELFF_AMREX_HYDRO_PREFIX "${CMAKE_INSTALL_PREFIX}" CACHE PATH "AMReX-Hydro install prefix")
        add_subdirectory(amrex-hydro)
    endif()
endif()

if (ELFF_BUILD_TESTS)
    if(ELFF_BUILD_THIRD_PARTY_GTEST)
        set(ELFF_GTEST_PREFIX "${CMAKE_INSTALL_PREFIX}" CACHE PATH "Gtest install prefix")
        add_subdirectory(gtest)
    endif()
endif()
