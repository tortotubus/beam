
set(CMAKE_VERBOSE_MAKEFILE ON)

# Enable C language
enable_language(C)
 
elff_basilisk_mpi_add_executable("${CMAKE_CURRENT_SOURCE_DIR}/mpi_huang.c")
elff_basilisk_mpi_add_executable("${CMAKE_CURRENT_SOURCE_DIR}/mpi_uhllman_cylinder.c")
elff_basilisk_add_executable("${CMAKE_CURRENT_SOURCE_DIR}/huang.c")
elff_basilisk_add_executable("${CMAKE_CURRENT_SOURCE_DIR}/uhllman_cylinder.c")