enable_testing()

file(GLOB BENCHMARK_TEST_SOURCES benchmark/*.cpp benchmark/*.hpp)
add_executable(elff_benchmark_tests ${BENCHMARK_TEST_SOURCES})
target_include_directories(elff_benchmark_tests PRIVATE ${CMAKE_SOURCE_DIR})
target_link_libraries(elff_benchmark_tests PRIVATE ELFF GTest::gtest_main)

if(ELFF_BUILD_THIRD_PARTY_GTEST)
    add_test(
        NAME EulerBeamDynamicInextensibleMoMTest
        COMMAND elff_benchmark_tests --gtest_filter=EulerBeamDynamicInextensibleMoMTest.*
    )
    add_test(
        NAME EulerBeamDynamicInextensibleMoMSparseTest
        COMMAND elff_benchmark_tests --gtest_filter=EulerBeamDynamicInextensibleMoMSparseTest.*
    )
    add_test(
        NAME EulerBeamStaticInextensibleMoMTest
        COMMAND elff_benchmark_tests --gtest_filter=EulerBeamStaticInextensibleMoMTest.*
    )
    add_test(
        NAME EulerBeamStaticInextensibleMoMSparseTest
        COMMAND elff_benchmark_tests --gtest_filter=EulerBeamStaticInextensibleMoMSparseTest.*
    )
else()
    gtest_discover_tests(elff_benchmark_tests)
endif()