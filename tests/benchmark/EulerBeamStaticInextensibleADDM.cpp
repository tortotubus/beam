#include <elff/models/beam/EulerBeamStaticInextensibleADDM.hpp>
#include <gtest/gtest.h>

#include "EulerBeamStaticInextensibleReferences.hpp"

namespace ELFF {
using namespace IO::CXX;
using namespace Models;

TEST(EulerBeamStaticInextensibleADDMTest, BisshoppAndDrucker)
{

  real_t length = 1., EI = 1., area = 1., r_pentalty = 1e3;
  size_t nodes = 40;

  real_t tip_force_y = -1;

  double comparison_tol = 1e-3;

  EulerBeam::EulerBeamBCs boundary_conditions = {
    .end = { EulerBeam::left, EulerBeam::right },
    .type = { EulerBeam::clamped_bc, EulerBeam::point_force_bc },
    .vals = { { .position = { 0, 0, 0 }, .slope = { 1, 0, 0 } },
              { .force = { 0, tip_force_y, 0 } } }
  };

  EulerBeamStaticInextensibleADDM static_beam(
    length, EI, nodes, boundary_conditions, r_pentalty);
  static_beam.solve();

  EulerBeamMesh& mesh = static_beam.get_mesh();
  auto centerline = mesh.get_centerline();
  std::array<real_t, 3> tip = centerline[nodes - 1];

  BisshoppAndDrucker1945Result res =
    BisshoppAndDrucker1945(length, EI, -tip_force_y);

  EXPECT_NEAR(std::abs(length - tip[0]), res.A, comparison_tol);
  EXPECT_NEAR(std::abs(tip[1]), res.delta, comparison_tol);
};

}
