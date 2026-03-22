enable_testing()

# file(GLOB_RECURSE UNIT_TEST_SOURCES unit/*.cpp)
# add_executable(elff_tests ${UNIT_TEST_SOURCES})
# target_include_directories(elff_tests PRIVATE ${CMAKE_SOURCE_DIR})
# target_link_libraries(elff_tests PRIVATE elff GTest::gtest_main)
# gtest_discover_tests(elff_tests)


file(GLOB BENCHMARK_TEST_SOURCES benchmark/*.cpp benchmark/*.hpp)
add_executable(elff_benchmark_tests ${BENCHMARK_TEST_SOURCES})
target_include_directories(elff_benchmark_tests PRIVATE ${CMAKE_SOURCE_DIR})
target_link_libraries(elff_benchmark_tests PRIVATE elff GTest::gtest_main)
gtest_discover_tests(elff_benchmark_tests)
