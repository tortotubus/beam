# tests/CMakeLists.txt
enable_testing()

# Point it at your test sources
file(GLOB_RECURSE UNIT_TEST_SOURCES unit/*.cpp)
add_executable(beam_unit_tests ${UNIT_TEST_SOURCES})
target_include_directories(beam_unit_tests PRIVATE ${CMAKE_SOURCE_DIR})
target_link_libraries(beam_unit_tests PRIVATE beam GTest::gtest_main)
gtest_discover_tests(beam_unit_tests)


file(GLOB BENCHMARK_TEST_SOURCES benchmark/*.cpp benchmark/*.hpp)
add_executable(beam_benchmark_tests ${BENCHMARK_TEST_SOURCES})
target_include_directories(beam_benchmark_tests PRIVATE ${CMAKE_SOURCE_DIR})
target_link_libraries(beam_benchmark_tests PRIVATE beam GTest::gtest_main)
gtest_discover_tests(beam_benchmark_tests)

# …and similarly for performance/ if desired
