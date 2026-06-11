
set(CMAKE_VERBOSE_MAKEFILE ON)

# Enable C language
enable_language(C)

elff_basilisk_add_executable("${CMAKE_CURRENT_SOURCE_DIR}/adapt_wavelet_ibm_test.c")
 
elff_basilisk_add_executable("${CMAKE_CURRENT_SOURCE_DIR}/huang_quadtree.c")
elff_basilisk_mpi_add_executable("${CMAKE_CURRENT_SOURCE_DIR}/huang_quadtree_mpi.c")

elff_basilisk_add_executable("${CMAKE_CURRENT_SOURCE_DIR}/rbc_shear.c")

elff_basilisk_add_executable("${CMAKE_CURRENT_SOURCE_DIR}/static_cylinder.c")
elff_basilisk_mpi_add_executable("${CMAKE_CURRENT_SOURCE_DIR}/static_cylinder_mpi.c")

elff_basilisk_add_executable("${CMAKE_CURRENT_SOURCE_DIR}/static_cylinder_3d_embed.c")
# elff_basilisk_mpi_add_executable("${CMAKE_CURRENT_SOURCE_DIR}/static_cylinder_mpi.c")

elff_basilisk_add_executable("${CMAKE_CURRENT_SOURCE_DIR}/static_sphere.c")
elff_basilisk_mpi_add_executable("${CMAKE_CURRENT_SOURCE_DIR}/static_sphere_mpi.c")

# elff_basilisk_add_executable("${CMAKE_CURRENT_SOURCE_DIR}/static_fibre.c")
# elff_basilisk_mpi_add_executable("${CMAKE_CURRENT_SOURCE_DIR}/static_fibre_mpi.c")

elff_basilisk_add_executable("${CMAKE_CURRENT_SOURCE_DIR}/single_settling_fibre.c")
elff_basilisk_mpi_add_executable("${CMAKE_CURRENT_SOURCE_DIR}/single_settling_fibre_mpi.c")
