enable_testing()
include(GoogleTest)

file(GLOB BENCHMARK_TEST_SOURCES CONFIGURE_DEPENDS benchmark/*.cpp benchmark/*.hpp)
add_executable(elff_benchmark_tests ${BENCHMARK_TEST_SOURCES})
target_include_directories(elff_benchmark_tests PRIVATE ${CMAKE_SOURCE_DIR})
target_link_libraries(elff_benchmark_tests PRIVATE ELFF GTest::gtest_main)

if(ELFF_BUILD_THIRD_PARTY_GTEST)
    add_test(
        NAME BendingRodSegmentsImplicitBenchmark
        COMMAND elff_benchmark_tests --gtest_filter=BendingRodSegmentsImplicitBenchmark.*
    )
    add_test(
        NAME BendingRodSegmentsIMEXBenchmark
        COMMAND elff_benchmark_tests --gtest_filter=BendingRodSegmentsIMEXBenchmark.*
    )
    add_test(
        NAME BendingRodSegmentsBenchmark
        COMMAND elff_benchmark_tests --gtest_filter=BendingRodSegmentsBenchmark.*
    )
    add_test(
        NAME EulerBeamDynamicMoMSparseTest
        COMMAND elff_benchmark_tests --gtest_filter=EulerBeamDynamicMoMSparseTest.*
    )
    add_test(
        NAME EulerBeamDynamicIndex1Test
        COMMAND elff_benchmark_tests --gtest_filter=EulerBeamDynamicIndex1Test.*
    )
    add_test(
        NAME EulerBeamDynamicGGLTest
        COMMAND elff_benchmark_tests --gtest_filter=EulerBeamDynamicGGLTest.*
    )
    add_test(
        NAME EulerBeamGGL2Test
        COMMAND elff_benchmark_tests --gtest_filter=EulerBeamGGL2Test.*
    )
    add_test(
        NAME EulerBeamHuangTest
        COMMAND elff_benchmark_tests --gtest_filter=EulerBeamHuangTest.*
    )
    add_test(
        NAME EulerBeamMoMTest
        COMMAND elff_benchmark_tests --gtest_filter=EulerBeamMoMTest.*
    )
    add_test(
        NAME EulerBeamStaticMoMSparseTest
        COMMAND elff_benchmark_tests --gtest_filter=EulerBeamStaticMoMSparseTest.*
    )
else()
    gtest_discover_tests(elff_benchmark_tests)
endif()

file(GLOB_RECURSE UNIT_TEST_SOURCES CONFIGURE_DEPENDS unit/*.cpp unit/*.hpp)
add_executable(elff_unit_tests ${UNIT_TEST_SOURCES})
target_include_directories(elff_unit_tests PRIVATE
    ${CMAKE_SOURCE_DIR} "${CMAKE_BINARY_DIR}/generated")
target_link_libraries(elff_unit_tests PRIVATE ELFF GTest::gtest_main)
target_compile_definitions(
    elff_unit_tests
    PRIVATE
        ELFF_CASE_DTD_PATH="${ELFF_V1_DTD_SOURCE}"
        ELFF_FLAPPING_FIBRE_CASE_PATH="${CMAKE_SOURCE_DIR}/basilisk/examples/cases/flapping_fibre.xml"
)

if(ELFF_BUILD_THIRD_PARTY_GTEST)
    add_test(NAME XMLNodeTest COMMAND elff_unit_tests --gtest_filter=XMLNodeTest.*)
    add_test(NAME EulerBeamHuangXMLTest COMMAND elff_unit_tests --gtest_filter=EulerBeamHuangXMLTest.*)
    add_test(
        NAME BasiliskCaseReaderTest
        COMMAND elff_unit_tests --gtest_filter=BasiliskCaseReaderTest.*
    )
    add_test(
        NAME ELFFCaseDTDTest
        COMMAND elff_unit_tests --gtest_filter=ELFFCaseDTDTest.*
    )
    add_test(
        NAME XMLCaseDocumentTest
        COMMAND elff_unit_tests --gtest_filter=XMLCaseDocument*.*
    )
    add_test(
        NAME BendingRodSegmentsImplicitTest
        COMMAND elff_unit_tests --gtest_filter=BendingRodSegmentsImplicitTest.*
    )
    add_test(
        NAME BendingRodSegmentsIMEXTest
        COMMAND elff_unit_tests --gtest_filter=BendingRodSegmentsIMEXTest.*
    )
    add_test(
        NAME BendingRodSegmentsTest
        COMMAND elff_unit_tests --gtest_filter=BendingRodSegmentsTest.*
    )
    add_test(
        NAME vtkHDFPolyDataStaticWriterTest
        COMMAND elff_unit_tests --gtest_filter=vtkHDFPolyDataStaticWriterTest.*
    )
    add_test(
        NAME EulerBeamDynamicGGLUnitTest
        COMMAND elff_unit_tests --gtest_filter=EulerBeamDynamicGGLUnitTest.*
    )
else()
    gtest_discover_tests(elff_unit_tests)
endif()

if(TARGET periodic_fibre_serial)
    add_test(
        NAME BasiliskPeriodicFibreSerial
        COMMAND $<TARGET_FILE:periodic_fibre_serial>
    )
    set_tests_properties(BasiliskPeriodicFibreSerial PROPERTIES TIMEOUT 60)
endif()

if(TARGET flapping_fibre_serial)
    add_test(
        NAME BasiliskFlappingFibreSerial
        COMMAND $<TARGET_FILE:flapping_fibre_serial>
                "${CMAKE_SOURCE_DIR}/tests/data/flapping_fibre_short.xml"
    )
    set_tests_properties(
        BasiliskFlappingFibreSerial
        PROPERTIES
            TIMEOUT 60
            WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    )
endif()

if(TARGET periodic_fibre_mpi AND MPIEXEC_EXECUTABLE)
    add_test(
        NAME BasiliskPeriodicFibreMPI
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                ${MPIEXEC_PREFLAGS} $<TARGET_FILE:periodic_fibre_mpi>
                ${MPIEXEC_POSTFLAGS}
    )
    set_tests_properties(
        BasiliskPeriodicFibreMPI PROPERTIES PROCESSORS 2 TIMEOUT 60
    )
endif()

if(TARGET flapping_fibre_mpi AND MPIEXEC_EXECUTABLE)
    add_test(
        NAME BasiliskFlappingFibreMPI
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                ${MPIEXEC_PREFLAGS} $<TARGET_FILE:flapping_fibre_mpi>
                "${CMAKE_SOURCE_DIR}/tests/data/flapping_fibre_short.xml"
                ${MPIEXEC_POSTFLAGS}
    )
    set_tests_properties(
        BasiliskFlappingFibreMPI
        PROPERTIES
            PROCESSORS 2
            TIMEOUT 60
            WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    )
endif()

find_package(Python3 3.7 REQUIRED COMPONENTS Interpreter)
add_test(NAME ModelDTDGeneration
    COMMAND "${Python3_EXECUTABLE}"
        "${CMAKE_CURRENT_SOURCE_DIR}/python/test_bundle_dtd.py")
