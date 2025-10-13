

#include <gtest/gtest.h>
#include "models/beam/EulerBeamStaticInextensibleMoMSparse.hpp"

using beam::EulerBeamStaticInextensibleMoMSparse;
using beam::EulerBeam;
using beam::EulerBeamBCs;
using beam::real_t;
using beam::left;
using beam::right;
using beam::clamped_bc;
using beam::free_bc;

TEST(EulerBeamStaticInextensibleMoMSparseTest, DumpLoadAndStiffness) {

  EulerBeamBCs boundary_conditions = {
    .end = {left, right},
    .type = {clamped_bc, free_bc},
    .vals = {{
      .position = {0,0,0},
      .slope = {1,0,0}
    }, {
      .position = {1,0,0},
      .slope = {1,0,0}
    }}
  };

  real_t length = 1., EI = 1., load = -1., area = 1., r_penalty = 1e4;
  size_t nodes = 3;

  EulerBeamStaticInextensibleMoMSparse static_beam(length, EI, load, area, nodes, boundary_conditions, r_penalty);

  //static_beam.update_lambda();
  static_beam.apply_initial_condition();
  static_beam.assemble_system();

  static_beam.apply_boundary_conditions();

  auto J = static_beam.jacobian;
  auto u = static_beam.u;
  auto R = static_beam.residual;
  
  GTEST_LOG_(INFO) << "CTEST_FULL_OUTPUT";
  GTEST_LOG_(INFO) << "J =\n" << J << "\n";
  GTEST_LOG_(INFO) << "u =\n" << u << "\n";
  GTEST_LOG_(INFO) << "R =\n" << R << "\n";
  
  // Optionally, add some real checks:
  // EXPECT_NEAR(F(0), expected_value, 1e-12);
  // EXPECT_NEAR(K(0,0), expected_K00, 1e-12);
};


TEST(EulerBeamStaticInextensibleMoMSparseTest, SolveUniformLoadAndPlot) {

  real_t length = 1., EI = 1., area = 1., r_pentalty = 1e4, load = -10;
  size_t nodes = 100;

  EulerBeamBCs boundary_conditions = {
    .end = {left, right},
    .type = {clamped_bc, free_bc},
    .vals = {{
      .position = {0,0,0},
      .slope = {1,0,0}
    }, {
      .position = {1,0,0},
      .slope = {1,0,0}
    }}
  };

  EulerBeamStaticInextensibleMoMSparse static_beam(length, EI, load, area, nodes, boundary_conditions, r_pentalty);

  GTEST_LOG_(INFO) << "CTEST_FULL_OUTPUT";
  static_beam.solve();
  static_beam.plot("Static Inextensible Euler Beam");
  
  // auto u = static_beam.get_solution();
  // GTEST_LOG_(INFO) << "u =\n" << u << "\n";

  // Optionally, add some real checks:
  // EXPECT_NEAR(F(0), expected_value, 1e-12);
  // EXPECT_NEAR(K(0,0), expected_K00, 1e-12);
};