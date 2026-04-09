
set(CMAKE_VERBOSE_MAKEFILE ON)

# Enable C language
enable_language(C)
 
elff_basilisk_mpi_add_executable("${CMAKE_CURRENT_SOURCE_DIR}/huang_quadtree_mpi.c")
# elff_basilisk_mpi_add_executable("${CMAKE_CURRENT_SOURCE_DIR}/huang_two_fibres.c")
# elff_basilisk_mpi_add_executable("${CMAKE_CURRENT_SOURCE_DIR}/uhllman_cylinder.c")

elff_basilisk_add_executable("${CMAKE_CURRENT_SOURCE_DIR}/huang_quadtree.c")
elff_basilisk_add_executable("${CMAKE_CURRENT_SOURCE_DIR}/huang_quadtree_two_fibres.c")
elff_basilisk_add_executable("${CMAKE_CURRENT_SOURCE_DIR}/uhllman_cylinder.c")