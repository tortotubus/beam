#include <elff/models/beam/EulerBeamStaticInextensibleADDM.hpp>
#include <gtest/gtest.h>

#include "EulerBeamStaticInextensibleReferences.hpp"

namespace ELFF {
using namespace IO::CXX;
using namespace Models;

TEST(EulerBeamStaticInextensibleADDMTest, BisshoppAndDrucker)
{

  real_t length = 1., EI = 1., area = 1., r_pentalty = 1e2;
  size_t nodes = 40;

  real_t tip_force_y = -1;

  double comparison_tol = 5e-7;

  EulerBeam::EulerBeamBCs boundary_conditions = {
    .end = { EulerBeam::left, EulerBeam::right },
    .type = { EulerBeam::clamped_bc, EulerBeam::point_force_bc },
    .vals = { { .position = { 0, 0, 0 }, .slope = { 1, 0, 0 } },
              { .force = { 0, tip_force_y, 0 } } }
  };

  EulerBeamStaticInextensibleADDM static_beam(
    length, EI, nodes, boundary_conditions, r_pentalty);
  static_beam.solve();
  static_beam.get_mesh().plot_gnuplot("Bisshopp and Drucker ADDM");

  EulerBeamMesh& mesh = static_beam.get_mesh();
  auto centerline = mesh.get_centerline();
  std::array<real_t, 3> tip = centerline[nodes - 1];

  BisshoppAndDrucker1945Result res =
    BisshoppAndDrucker1945(length, EI, -tip_force_y);

  EXPECT_NEAR(std::abs(length - tip[0]), res.A, comparison_tol);
  EXPECT_NEAR(std::abs(tip[1]), res.delta, comparison_tol);
};

TEST(EulerBeamStaticInextensibleADDMTest, GlowinskiIC)
{
  GTEST_LOG_(INFO) << "CTEST_FULL_OUTPUT";
  real_t length = 32.6, EI = 700., mu = 7.67, r_penalty = 1e2;
  std::array<real_t, 3> load = { 0, -9.81 * mu, 0 };

  real_t dt = 1e-2;
  real_t t = 0;
  real_t tf = 10;
  size_t Nt = size_t(ceil(tf / dt));

  real_t dt_save = 0.1;
  size_t Nt_save = size_t(ceil(dt_save / dt));

  size_t nodes = 61;

  EulerBeam::EulerBeamBCs boundary_conditions = {
    .end = { EulerBeam::left, EulerBeam::right },
    .type = { EulerBeam::simple_bc, EulerBeam::simple_bc },
    .vals = { {
                .position = { 0, 0, 0 },
              },
              {
                .position = { 20, 0, 0 },
              } }
  };

  EulerBeamStaticInextensibleADDM static_beam(
    length, EI, nodes, boundary_conditions, r_penalty);
  static_beam.apply_initial_condition();

  ELFF_LOG("Static Solve:");
  static_beam.solve(load);
  static_beam.get_mesh().plot_gnuplot("Glowinski IC ADDM");
};
}
