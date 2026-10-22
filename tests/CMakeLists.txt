enable_testing()
include(GoogleTest)

file(GLOB BENCHMARK_TEST_SOURCES CONFIGURE_DEPENDS benchmark/*.cpp benchmark/*.hpp)
add_executable(elff_benchmark_tests ${BENCHMARK_TEST_SOURCES})
target_include_directories(elff_benchmark_tests PRIVATE ${CMAKE_SOURCE_DIR})
target_link_libraries(elff_benchmark_tests PRIVATE ELFF GTest::gtest_main) 

gtest_discover_tests(elff_benchmark_tests
  TEST_PREFIX "benchmark."
  DISCOVERY_MODE PRE_TEST
)

