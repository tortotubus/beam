# tests/CMakeLists.txt
enable_testing()

# Point it at your test sources
file(GLOB_RECURSE UNIT_TEST_SOURCES unit/*.cpp)
add_executable(elff_unit_tests ${UNIT_TEST_SOURCES})
target_include_directories(elff_unit_tests PRIVATE ${CMAKE_SOURCE_DIR})
target_include_directories(elff_unit_tests PRIVATE ${CMAKE_SOURCE_DIR}/basilisk)
target_link_libraries(elff_unit_tests PRIVATE elff GTest::gtest_main)
gtest_discover_tests(elff_unit_tests)


file(GLOB BENCHMARK_TEST_SOURCES benchmark/*.cpp benchmark/*.hpp)
add_executable(elff_benchmark_tests ${BENCHMARK_TEST_SOURCES})
target_include_directories(elff_benchmark_tests PRIVATE ${CMAKE_SOURCE_DIR})
target_link_libraries(elff_benchmark_tests PRIVATE elff GTest::gtest_main)
gtest_discover_tests(elff_benchmark_tests)

# …and similarly for performance/ if desired
