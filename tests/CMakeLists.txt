enable_testing()
include(GoogleTest)

file(GLOB BENCHMARK_TEST_SOURCES CONFIGURE_DEPENDS benchmark/*.cpp benchmark/*.hpp)
add_executable(elff_benchmark_tests ${BENCHMARK_TEST_SOURCES})
target_include_directories(elff_benchmark_tests PRIVATE ${CMAKE_SOURCE_DIR})
target_link_libraries(elff_benchmark_tests PRIVATE ELFF GTest::gtest_main)

if(ELFF_BUILD_THIRD_PARTY_GTEST)
    add_test(
        NAME EulerBeamDynamicInextensibleMoMSparseTest
        COMMAND elff_benchmark_tests --gtest_filter=EulerBeamDynamicInextensibleMoMSparseTest.*
    )
    add_test(
        NAME EulerBeamDynamicInextensibleIndex1Test
        COMMAND elff_benchmark_tests --gtest_filter=EulerBeamDynamicInextensibleIndex1Test.*
    )
    add_test(
        NAME EulerBeamDynamicInextensibleGGLTest
        COMMAND elff_benchmark_tests --gtest_filter=EulerBeamDynamicInextensibleGGLTest.*
    )
    add_test(
        NAME EulerBeamInextensibleGGL2Test
        COMMAND elff_benchmark_tests --gtest_filter=EulerBeamInextensibleGGL2Test.*
    )
    add_test(
        NAME EulerBeamHuangTest
        COMMAND elff_benchmark_tests --gtest_filter=EulerBeamHuangTest.*
    )
    add_test(
        NAME EulerBeamInextensibleMoMTest
        COMMAND elff_benchmark_tests --gtest_filter=EulerBeamInextensibleMoMTest.*
    )
    add_test(
        NAME EulerBeamStaticInextensibleMoMSparseTest
        COMMAND elff_benchmark_tests --gtest_filter=EulerBeamStaticInextensibleMoMSparseTest.*
    )
else()
    gtest_discover_tests(elff_benchmark_tests)
endif()

file(GLOB_RECURSE UNIT_TEST_SOURCES CONFIGURE_DEPENDS unit/*.cpp unit/*.hpp)
add_executable(elff_unit_tests ${UNIT_TEST_SOURCES})
target_include_directories(elff_unit_tests PRIVATE ${CMAKE_SOURCE_DIR})
target_link_libraries(elff_unit_tests PRIVATE ELFF GTest::gtest_main)

if(ELFF_BUILD_THIRD_PARTY_GTEST)
    add_test(
        NAME vtkHDFPolyDataStaticWriterTest
        COMMAND elff_unit_tests --gtest_filter=vtkHDFPolyDataStaticWriterTest.*
    )
    add_test(
        NAME EulerBeamDynamicInextensibleGGLUnitTest
        COMMAND elff_unit_tests --gtest_filter=EulerBeamDynamicInextensibleGGLUnitTest.*
    )
else()
    gtest_discover_tests(elff_unit_tests)
endif()
