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

# Focused runtime composition and migrated public API regressions.
add_executable(elff_runtime_tests
  unit/runtime/Composition.cpp
  unit/models/ModelRuntime.cpp
  unit/models/ModelRuntimeXML.cpp
  unit/models/BendingRodImplicitXML.cpp
  unit/models/beam/IBEulerBeamCAPI.cpp
  unit/models/capsule/IBCapsuleCAPI.cpp
  unit/models/rigidbodies/IBRigidBodyCAPI.cpp
  unit/models/ibm/IBSpringCAPI.cpp)
target_link_libraries(elff_runtime_tests PRIVATE ELFF GTest::gtest_main)
target_include_directories(elff_runtime_tests PRIVATE "${CMAKE_BINARY_DIR}/generated")
target_compile_definitions(elff_runtime_tests PRIVATE
  ELFF_LEGACY_RUNTIME_CHECKPOINT="${CMAKE_CURRENT_SOURCE_DIR}/data/model_runtime_v2_legacy.bin")
add_test(NAME RuntimeCompositionTest COMMAND elff_runtime_tests)

configure_file(data/runtime_composition.xml runtime_composition.xml COPYONLY)

if(ELFF_USE_MPI AND ELFF_USE_LIBXML2)
  add_executable(elff_runtime_communicator_mpi mpi/ModelRuntimeCommunicator.cpp)
  add_executable(elff_runtime_implicit_rod_mpi mpi/ModelRuntimeImplicitRod.cpp)
  target_link_libraries(elff_runtime_communicator_mpi PRIVATE ELFF)
  target_link_libraries(elff_runtime_implicit_rod_mpi PRIVATE ELFF)
  add_test(NAME ModelRuntimeCommunicatorMPI COMMAND ${MPIEXEC_EXECUTABLE}
    ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS} $<TARGET_FILE:elff_runtime_communicator_mpi>
    ${MPIEXEC_POSTFLAGS} "${CMAKE_CURRENT_BINARY_DIR}/runtime_composition.xml")
  add_test(NAME ModelRuntimeImplicitRodMPI COMMAND ${MPIEXEC_EXECUTABLE}
    ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS} $<TARGET_FILE:elff_runtime_implicit_rod_mpi>
    ${MPIEXEC_POSTFLAGS})
  set_tests_properties(ModelRuntimeCommunicatorMPI ModelRuntimeImplicitRodMPI
    PROPERTIES PROCESSORS 2 TIMEOUT 60)
endif()

if(ELFF_USE_LIBXML2 AND TARGET flapping_fibre_serial)
  file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/basilisk_serial")
  configure_file(data/runtime_composition.xml basilisk_serial/case.xml COPYONLY)
  add_test(NAME ModelRuntimeBasiliskSerial COMMAND flapping_fibre_serial case.xml)
  set_tests_properties(ModelRuntimeBasiliskSerial PROPERTIES
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/basilisk_serial" TIMEOUT 60)
endif()
if(ELFF_USE_MPI AND ELFF_USE_LIBXML2 AND TARGET flapping_fibre_mpi)
  file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/basilisk_mpi")
  configure_file(data/runtime_composition.xml basilisk_mpi/case.xml COPYONLY)
  add_test(NAME ModelRuntimeBasiliskMPI COMMAND ${MPIEXEC_EXECUTABLE}
    ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS} $<TARGET_FILE:flapping_fibre_mpi>
    ${MPIEXEC_POSTFLAGS} case.xml)
  set_tests_properties(ModelRuntimeBasiliskMPI PROPERTIES
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/basilisk_mpi" PROCESSORS 2 TIMEOUT 60)
endif()
