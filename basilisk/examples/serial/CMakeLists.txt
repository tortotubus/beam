
set(CMAKE_VERBOSE_MAKEFILE ON)

# Enable C language
enable_language(C)

elff_basilisk_add_executable("${CMAKE_CURRENT_SOURCE_DIR}/taylor-green-basilisk.c")
elff_basilisk_add_executable("${CMAKE_CURRENT_SOURCE_DIR}/taylor-green.c")
