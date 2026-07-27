
set(CMAKE_VERBOSE_MAKEFILE ON)

# Enable C language
enable_language(C)

elff_basilisk_add_executable("${CMAKE_CURRENT_SOURCE_DIR}/adapt_wavelet_ibm_test.c")

elff_basilisk_add_executable(
  "${CMAKE_CURRENT_SOURCE_DIR}/decaying_taylor_vortex_fibres.c"
  -disable-dimensions)
elff_basilisk_mpi_add_executable(
  "${CMAKE_CURRENT_SOURCE_DIR}/decaying_taylor_vortex_fibres_mpi.c"
  -disable-dimensions)
 
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

elff_basilisk_add_executable("${CMAKE_CURRENT_SOURCE_DIR}/single_settling_fibre_admm.c")
elff_basilisk_mpi_add_executable("${CMAKE_CURRENT_SOURCE_DIR}/single_settling_fibre_admm_mpi.c")

elff_basilisk_add_executable("${CMAKE_CURRENT_SOURCE_DIR}/single_settling_fibre_huang.c")
elff_basilisk_mpi_add_executable("${CMAKE_CURRENT_SOURCE_DIR}/single_settling_fibre_huang_mpi.c")

elff_basilisk_add_executable("${CMAKE_CURRENT_SOURCE_DIR}/many_settling_fibres_huang.c")
elff_basilisk_mpi_add_executable("${CMAKE_CURRENT_SOURCE_DIR}/many_settling_fibres_huang_mpi.c")

elff_basilisk_add_executable("${CMAKE_CURRENT_SOURCE_DIR}/single_settling_fibre_ggl.c")
elff_basilisk_mpi_add_executable("${CMAKE_CURRENT_SOURCE_DIR}/single_settling_fibre_ggl_mpi.c")

elff_basilisk_add_executable("${CMAKE_CURRENT_SOURCE_DIR}/quadtree_moving_mesh.c")
elff_basilisk_mpi_add_executable("${CMAKE_CURRENT_SOURCE_DIR}/quadtree_moving_mesh_mpi.c")
