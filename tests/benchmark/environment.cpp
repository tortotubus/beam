// test_env.cpp
#include <gtest/gtest.h>
#include <Eigen/Core>
#include <omp.h>

struct EigenEnv : ::testing::Environment {
  void SetUp() override {
    // Size both OpenMP and Eigen once for the whole test process
    const int n = omp_get_max_threads();
    omp_set_num_threads(n);
    Eigen::setNbThreads(n);
  }
};

// Static registration runs before gtest_main's RUN_ALL_TESTS()
::testing::Environment* const eigen_env =
    ::testing::AddGlobalTestEnvironment(new EigenEnv);
