# tests/CMakeLists.txt
enable_testing()

# Find or FetchContent GoogleTest here (as you did in top‐level)
# add_subdirectory(...) if you vendor GTest under tests/

# Point it at your test sources
file(GLOB_RECURSE UNIT_TEST_SOURCES unit/*.cpp)
#file(GLOB INT_TEST_SOURCES integration/*.cpp)
#file(GLOB PERF_TEST_SOURCES performance/*.cpp)

add_executable(beam_unit_tests ${UNIT_TEST_SOURCES})
target_include_directories(beam_unit_tests PRIVATE ${CMAKE_SOURCE_DIR})
target_link_libraries(beam_unit_tests PRIVATE beam GTest::gtest_main)
gtest_discover_tests(beam_unit_tests)

#add_executable(beam_integration_tests ${INT_TEST_SOURCES})
#target_link_libraries(beam_integration_tests PRIVATE beam GTest::gtest_main)
#gtest_discover_tests(beam_integration_tests)

# …and similarly for performance/ if desired
