
set(CMAKE_VERBOSE_MAKEFILE ON)

# Enable C language
enable_language(C)
 
elff_basilisk_mpi_add_executable("${CMAKE_CURRENT_SOURCE_DIR}/mpi_karman.c")
elff_basilisk_mpi_add_executable("${CMAKE_CURRENT_SOURCE_DIR}/mpi_io.c")
elff_basilisk_mpi_add_executable("${CMAKE_CURRENT_SOURCE_DIR}/mpi_kernels.c")
elff_basilisk_mpi_add_executable("${CMAKE_CURRENT_SOURCE_DIR}/mpi_vtkhdfpolydata.c")
elff_basilisk_mpi_add_executable("${CMAKE_CURRENT_SOURCE_DIR}/mpi_uhllman_cylinder.c")