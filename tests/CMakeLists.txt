enable_testing()
include(GoogleTest)

file(GLOB BENCHMARK_TEST_SOURCES CONFIGURE_DEPENDS benchmark/*.cpp benchmark/*.hpp)
add_executable(elff_benchmark_tests ${BENCHMARK_TEST_SOURCES})
target_include_directories(elff_benchmark_tests PRIVATE ${CMAKE_SOURCE_DIR})
target_link_libraries(elff_benchmark_tests PRIVATE ELFF GTest::gtest_main)
# Contact geometry checks are independent of model and XML fixtures.
file(GLOB CONTACT_TEST_SOURCES CONFIGURE_DEPENDS unit/contact/*.cpp)
add_executable(elff_contact_tests ${CONTACT_TEST_SOURCES})
target_link_libraries(elff_contact_tests PRIVATE ELFF GTest::gtest_main)
add_test(NAME ContactGeometryTest COMMAND elff_contact_tests)
