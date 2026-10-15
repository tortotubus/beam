#include "elff/models/rigidbody/RigidBodyCircle.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <vector>

namespace ELFF {
namespace Models {
namespace {

using Vec3 = RigidBodyModel::Vec3;

TEST(RigidBodyCircleTest, PerimeterGeometryAndMeasures)
{
  const real_t R = 2.0;
  const size_t N = 240;
  RigidBodyCircle circle(R, N, 1.0);

  const auto& pts = circle.mesh().reference_points();
  const auto& nrm = circle.mesh().reference_normals();
  const auto& ds = circle.mesh().ds();

  ASSERT_EQ(pts.size(), N);
  ASSERT_EQ(nrm.size(), N);
  ASSERT_EQ(ds.size(), N);

  const real_t pi = std::acos(-1.);
  const real_t expected_ds = 2. * pi * R / static_cast<real_t>(N);
  const real_t expected_perimeter = 2. * pi * R;

  for (size_t i = 0; i < N; ++i) {
    const real_t pnorm =
      std::sqrt(pts[i][0] * pts[i][0] + pts[i][1] * pts[i][1] + pts[i][2] * pts[i][2]);
    const real_t nnorm =
      std::sqrt(nrm[i][0] * nrm[i][0] + nrm[i][1] * nrm[i][1] + nrm[i][2] * nrm[i][2]);

    EXPECT_NEAR(pnorm, R, 1e-10);
    EXPECT_NEAR(pts[i][2], 0.0, 1e-12);
    EXPECT_NEAR(nnorm, 1.0, 1e-10);
    EXPECT_NEAR(nrm[i][2], 0.0, 1e-12);
    EXPECT_NEAR(ds[i], expected_ds, 1e-12);
  }

  EXPECT_NEAR(circle.mesh().total_measure(), expected_perimeter, 1e-10);
}

TEST(RigidBodyCircleTest, DensityControlsMassInStepResponse)
{
  const real_t R = 1.0;
  const size_t N = 128;
  const real_t rho = 2.0;
  RigidBodyCircle circle(R, N, rho);

  std::vector<Vec3> traction(N, { 1.0, 0.0, 0.0 });
  const real_t dt = 0.1;
  circle.step(dt, traction);

  const real_t pi = std::acos(-1.);
  const real_t m = rho * pi * R * R;
  const real_t F = 2. * pi * R;
  const real_t v_expected = dt * (F / m);

  EXPECT_NEAR(circle.com_velocity()[0], v_expected, 1e-10);
  EXPECT_NEAR(circle.com_velocity()[1], 0.0, 1e-12);
  EXPECT_NEAR(circle.com_velocity()[2], 0.0, 1e-12);
}

} // namespace
} // namespace Models
} // namespace ELFF
