

if (ELFF_BUILD_THIRD_PARTY_HDF5)
    set(ELFF_HDF5_PREFIX "${CMAKE_INSTALL_PREFIX}" CACHE PATH "HDF5 install prefix")
    add_subdirectory(hdf5)
endif()

if (ELFF_BUILD_THIRD_PARTY_BASILISK)
    set(ELFF_BASILISK_PREFIX "${CMAKE_INSTALL_PREFIX}" CACHE PATH "Basilisk install prefix")
    add_subdirectory(basilisk)
endif()

if (ELFF_BUILD_THIRD_PARTY_EIGEN3)
    set(ELFF_EIGEN3_PREFIX "${CMAKE_INSTALL_PREFIX}" CACHE PATH "Eigen3 install prefix")
    add_subdirectory(eigen)
endif()
