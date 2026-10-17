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
target_compile_definitions(
    elff_unit_tests
    PRIVATE
        ELFF_CASE_DTD_PATH="${CMAKE_SOURCE_DIR}/io/DTD/ELFF/ELFF-v1.dtd"
        ELFF_BASILISK_DTD_PATH="${CMAKE_SOURCE_DIR}/basilisk/templates/flapping_fibre-v1.dtd"
        ELFF_FLAPPING_FIBRE_CASE_PATH="${CMAKE_SOURCE_DIR}/basilisk/examples/cases/flapping_fibre.xml"
)

if(ELFF_BUILD_THIRD_PARTY_GTEST)
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
        NAME EulerBeamDynamicInextensibleGGLUnitTest
        COMMAND elff_unit_tests --gtest_filter=EulerBeamDynamicInextensibleGGLUnitTest.*
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
