#include "elff/models/rigidbody/RigidBodyFibre.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <vector>

namespace ELFF {
namespace Models {
namespace {

using Vec3 = RigidBodyFibre::Vec3;
using Quat = std::array<real_t, 4>;

class TestRigidBodyFibre : public RigidBodyFibre
{
public:
  using RigidBodyFibre::RigidBodyFibre;

  void set_state(const Vec3& x,
                 const Vec3& v,
                 const Quat& q,
                 const Vec3& omega)
  {
    set_initial_state(x, v, q, omega);
  }
};

TEST(RigidBodyFibreTest, CenterlineGeometryAndMeasures)
{
  const real_t L = 2.0;
  const real_t d = 0.1;
  const size_t N = 5;
  RigidBodyFibre fibre(L, d, N, 3.0);

  const auto& points = fibre.mesh().reference_points();
  const auto& ds = fibre.mesh().ds();
  ASSERT_EQ(points.size(), N);
  ASSERT_EQ(ds.size(), N);

  EXPECT_NEAR(points.front()[0], -L / 2.0, 1e-12);
  EXPECT_NEAR(points.back()[0], L / 2.0, 1e-12);
  EXPECT_NEAR(points[N / 2][0], 0.0, 1e-12);
  for (size_t i = 0; i < N; ++i) {
    EXPECT_NEAR(points[i][1], 0.0, 1e-12);
    EXPECT_NEAR(points[i][2], 0.0, 1e-12);
  }

  const real_t h = L / static_cast<real_t>(N - 1);
  EXPECT_NEAR(ds.front(), 0.5 * h, 1e-12);
  EXPECT_NEAR(ds.back(), 0.5 * h, 1e-12);
  EXPECT_NEAR(ds[1], h, 1e-12);
  EXPECT_NEAR(fibre.mesh().total_measure(), L, 1e-12);
}

TEST(RigidBodyFibreTest, MassAndCylinderInertia)
{
  const real_t L = 3.0;
  const real_t d = 0.2;
  const real_t lambda = 4.0;
  RigidBodyFibre fibre(L, d, 7, lambda);

  const real_t mass = lambda * L;
  const real_t r = 0.5 * d;
  const real_t Ixx = 0.5 * mass * r * r;
  const real_t Iyy = mass * (3.0 * r * r + L * L) / 12.0;

  EXPECT_NEAR(fibre.total_mass(), mass, 1e-12);
  EXPECT_NEAR(fibre.inertia_body()[0][0], Ixx, 1e-12);
  EXPECT_NEAR(fibre.inertia_body()[1][1], Iyy, 1e-12);
  EXPECT_NEAR(fibre.inertia_body()[2][2], Iyy, 1e-12);
  EXPECT_NEAR(fibre.inertia_body()[0][1], 0.0, 1e-12);
}

TEST(RigidBodyFibreTest, TransformIdentityAndRotation)
{
  TestRigidBodyFibre fibre(2.0, 0.1, 3, 1.0);

  EXPECT_NEAR(fibre.mesh().world_points().front()[0], -1.0, 1e-12);
  EXPECT_NEAR(fibre.mesh().world_points().back()[0], 1.0, 1e-12);

  const real_t c = std::sqrt(0.5);
  fibre.set_state(
    { 0., 0., 0. }, { 0., 0., 0. }, { c, 0., 0., c }, { 0., 0., 0. });

  EXPECT_NEAR(fibre.mesh().world_points().front()[0], 0.0, 1e-12);
  EXPECT_NEAR(fibre.mesh().world_points().front()[1], -1.0, 1e-12);
  EXPECT_NEAR(fibre.mesh().world_points().back()[0], 0.0, 1e-12);
  EXPECT_NEAR(fibre.mesh().world_points().back()[1], 1.0, 1e-12);
}

TEST(RigidBodyFibreDynamicsTest, ZeroTractionPreservesMomenta)
{
  TestRigidBodyFibre fibre(2.0, 0.1, 5, 2.0);
  fibre.set_state({ 1., 2., 3. },
                  { 0.5, -0.25, 0.75 },
                  { 1., 0., 0., 0. },
                  { 0.1, 0., 0. });

  const Vec3 v0 = fibre.com_velocity();
  const Vec3 L0 = fibre.angular_momentum_body();

  fibre.step(0.1, std::vector<Vec3>(fibre.nodes(), { 0., 0., 0. }));

  for (size_t k = 0; k < 3; ++k) {
    EXPECT_NEAR(fibre.com_velocity()[k], v0[k], 1e-12);
    EXPECT_NEAR(fibre.angular_momentum_body()[k], L0[k], 1e-12);
  }
}

TEST(RigidBodyFibreDynamicsTest, UniformForceDensityAcceleratesCom)
{
  const real_t L = 2.0;
  const real_t lambda = 4.0;
  RigidBodyFibre fibre(L, 0.1, 5, lambda);

  fibre.step(0.2, std::vector<Vec3>(fibre.nodes(), { 1., 0., 0. }));

  const real_t expected_vx = 0.2 * L / (lambda * L);
  EXPECT_NEAR(fibre.com_velocity()[0], expected_vx, 1e-12);
  EXPECT_NEAR(fibre.com_velocity()[1], 0.0, 1e-12);
  EXPECT_NEAR(fibre.com_velocity()[2], 0.0, 1e-12);
}

TEST(RigidBodyFibreDynamicsTest, AntisymmetricForceDensityProducesTorque)
{
  const real_t L = 2.0;
  const real_t dt = 0.05;
  RigidBodyFibre fibre(L, 0.2, 2, 3.0);

  std::vector<Vec3> traction = { { 0., -1., 0. }, { 0., 1., 0. } };
  fibre.step(dt, traction);

  const real_t torque_z = L * L / 2.0;
  const real_t expected_omega_z = dt * torque_z / fibre.inertia_body()[2][2];
  EXPECT_NEAR(fibre.angular_velocity_body()[0], 0.0, 1e-12);
  EXPECT_NEAR(fibre.angular_velocity_body()[1], 0.0, 1e-12);
  EXPECT_NEAR(fibre.angular_velocity_body()[2], expected_omega_z, 1e-12);
}

} // namespace
} // namespace Models
} // namespace ELFF
