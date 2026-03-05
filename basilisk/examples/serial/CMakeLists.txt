
set(CMAKE_VERBOSE_MAKEFILE ON)

# Enable C language
enable_language(C)
 
elff_basilisk_add_executable("${CMAKE_CURRENT_SOURCE_DIR}/taylor-green.c")
elff_basilisk_add_executable("${CMAKE_CURRENT_SOURCE_DIR}/taylor-green-centered.c")
elff_basilisk_add_executable("${CMAKE_CURRENT_SOURCE_DIR}/karman.c")
elff_basilisk_add_executable("${CMAKE_CURRENT_SOURCE_DIR}/kernels.c")
elff_basilisk_add_executable("${CMAKE_CURRENT_SOURCE_DIR}/uhllman_cylinder.c")
elff_basilisk_add_executable("${CMAKE_CURRENT_SOURCE_DIR}/vtkhdfpolydata.c")
