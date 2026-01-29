set(CMAKE_VERBOSE_MAKEFILE ON)

# Enable C language
enable_language(C)

find_package(PkgConfig REQUIRED)
pkg_check_modules(CRITERION REQUIRED criterion)

elff_basilisk_add_executable("${CMAKE_CURRENT_SOURCE_DIR}/basilisk_tests.c")
