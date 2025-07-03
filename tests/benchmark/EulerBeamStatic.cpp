#include <gtest/gtest.h>
#include "models/beam/EulerBeamStatic.hpp"

using beam::EulerBeam;
using beam::EulerBeamBCs;
using beam::EulerBeamMesh;
using beam::EulerBeamStatic;
using beam::real_t;
using beam::left;
using beam::right;
using beam::clamped_bc;
using beam::free_bc;

TEST(EulerBeamStaticTest, DumpLoadAndStiffness) {

  EulerBeamBCs boundary_conditions = {
    .end = {left, right},
    .type = {clamped_bc, free_bc},
    .vals = {{
      .position = {0,0,0},
      .slope = {0,0,0}
    }, {
      .position = {1,0,0},
      .slope = {0,0,0}
    }}
  };

  real_t length = 1., EI = 1., load = -1., area = 1.;
  size_t nodes = 3;

  EulerBeamStatic static_beam(length, EI, load, area, nodes, boundary_conditions);

  //static_beam.assemble_load();
  //static_beam.assemble_stiffness();
  //static_beam.apply_boundary_conditions();

  static_beam.solve();

  auto u = static_beam.u;
  auto F = static_beam.F_global;
  auto K = static_beam.K_global;

  GTEST_LOG_(INFO) << "F_global =\n" << F << "\n";
  GTEST_LOG_(INFO) << "K_global =\n" << K << "\n";
  GTEST_LOG_(INFO) << "u =\n" << u << "\n";
  
  // Optionally, add some real checks:
  // EXPECT_NEAR(F(0), expected_value, 1e-12);
  // EXPECT_NEAR(K(0,0), expected_K00, 1e-12);
};



TEST(EulerBeamStaticTest, SolveUniformLoadAndPlot) {

  real_t length = 1., EI = 1., load = -1., area = 1.;
  size_t nodes = 10;

  EulerBeamBCs boundary_conditions = {
    .end = {left, right},
    .type = {clamped_bc, free_bc},
    .vals = {{
      .position = {0,0,0},
      .slope = {0,0,0}
    }, {
      .position = {0,0,0},
      .slope = {0,0,0}
    }}
  };


  EulerBeamStatic static_beam(length, EI, load, area, nodes, boundary_conditions);

  static_beam.solve();

  GTEST_LOG_(INFO) << "F_global =\n" << static_beam.F_global << "\n";
  GTEST_LOG_(INFO) << "K_global =\n" << static_beam.K_global << "\n";
  GTEST_LOG_(INFO) << "u =\n" << static_beam.u << "\n";

  static_beam.plot();

  // Optionally, add some real checks:
  // EXPECT_NEAR(F(0), expected_value, 1e-12);
  // EXPECT_NEAR(K(0,0), expected_K00, 1e-12);
};


TEST(EulerBeamStaticTest, SolveNonUniformLoadAndPlot) {

  real_t length = 1., EI = 1., area = 1.;
  size_t nodes = 10;

  EulerBeamBCs boundary_conditions = {
    .end = {left, right},
    .type = {clamped_bc, free_bc},
    .vals = {{
      .position = {0,0,0},
      .slope = {0,0,0}
    }, {
      .position = {0,0,0},
      .slope = {0,0,0}
    }}
  };

  std::vector<std::array<real_t, 3>> load(nodes);

  for (size_t i = 0; i < nodes; i++) {
    load[i] = {1, -1, 0};
  }

  EulerBeamStatic static_beam(length, EI, load, area, nodes, boundary_conditions);

  static_beam.solve();

  GTEST_LOG_(INFO) << "F_global =\n" << static_beam.F_global << "\n";
  GTEST_LOG_(INFO) << "K_global =\n" << static_beam.K_global << "\n";
  GTEST_LOG_(INFO) << "u =\n" << static_beam.u << "\n";

  static_beam.plot();

  // Optionally, add some real checks:
  // EXPECT_NEAR(F(0), expected_value, 1e-12);
  // EXPECT_NEAR(K(0,0), expected_K00, 1e-12);
};

