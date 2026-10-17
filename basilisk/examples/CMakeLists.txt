
set(CMAKE_VERBOSE_MAKEFILE ON)

# Enable C language
enable_language(C)

elff_basilisk_add_executable("${CMAKE_CURRENT_SOURCE_DIR}/flapping_fibre_serial.c")
elff_basilisk_mpi_add_executable("${CMAKE_CURRENT_SOURCE_DIR}/flapping_fibre_mpi.c")

elff_basilisk_add_executable("${CMAKE_CURRENT_SOURCE_DIR}/single_settling_fibre_serial.c")
elff_basilisk_mpi_add_executable("${CMAKE_CURRENT_SOURCE_DIR}/single_settling_fibre_mpi.c")

elff_basilisk_add_executable("${CMAKE_CURRENT_SOURCE_DIR}/many_settling_fibres_serial.c")
elff_basilisk_mpi_add_executable("${CMAKE_CURRENT_SOURCE_DIR}/many_settling_fibres_mpi.c")
elff_basilisk_add_executable("${CMAKE_CURRENT_SOURCE_DIR}/isotropic_fibre_serial.c")
elff_basilisk_mpi_add_executable("${CMAKE_CURRENT_SOURCE_DIR}/isotropic_fibre_mpi.c")
