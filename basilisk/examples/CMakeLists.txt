
set(CMAKE_VERBOSE_MAKEFILE ON)

# Enable C language
enable_language(C)

elff_basilisk_add_executable("${CMAKE_CURRENT_SOURCE_DIR}/flapping_fibre_serial.c")
elff_basilisk_mpi_add_executable("${CMAKE_CURRENT_SOURCE_DIR}/flapping_fibre_mpi.c") 
elff_basilisk_add_executable("${CMAKE_CURRENT_SOURCE_DIR}/isotropic_fibre_serial.c")
elff_basilisk_mpi_add_executable("${CMAKE_CURRENT_SOURCE_DIR}/isotropic_fibre_mpi.c") 
elff_basilisk_add_executable("${CMAKE_CURRENT_SOURCE_DIR}/settling_fibre_serial.c")
elff_basilisk_mpi_add_executable("${CMAKE_CURRENT_SOURCE_DIR}/settling_fibre_mpi.c")

# Standalone communication benchmarks; not CTest tests.
# Synthetic ID-valued fields have no physical units; skip qcc dimensional interpretation.
elff_basilisk_mpi_add_executable("${CMAKE_CURRENT_SOURCE_DIR}/ib_scaling_quadtree_mpi.c" -disable-dimensions)
elff_basilisk_mpi_add_executable("${CMAKE_CURRENT_SOURCE_DIR}/ib_scaling_octree_mpi.c" -disable-dimensions)
