#include "elff/models/rigidbody/RigidBody2DFromReference.hpp"
#include "elff/general/error.hpp"

#include <cmath>
#include <gtest/gtest.h>

#include <vector>

namespace ELFF {
namespace Models {
namespace {

using Vec3 = RigidBodyModel::Vec3;

class TestRigidBody2DFromReference : public RigidBody2DFromReference
{
public:
  using RigidBody2DFromReference::RigidBody2DFromReference;

  void set_state(const Vec3& x, const Vec3& v, real_t theta, real_t omega_z)
  {
    set_initial_state(x, v, theta, omega_z);
  }
};

TEST(RigidBody2DMeshTest, RejectsInvalidDS)
{
  std::vector<Vec3> points = { { 0., 0., 0. }, { 1., 0., 0. } };
  std::vector<real_t> ds_bad = { 1., -1. };
  Vec3 cog = { 0., 0., 0. };

#ifdef ELFF_USE_EXCEPTIONS
  ELFF::set_error_action(ELFF::ELFF_ERROR_THROW);
  EXPECT_ANY_THROW({ (void)RigidBody2DFromReference(points, ds_bad, cog, 1., 1.); });
#else
  EXPECT_DEATH({ (void)RigidBody2DFromReference(points, ds_bad, cog, 1., 1.); },
               "DS entries must be positive");
#endif
}

TEST(RigidBody2DMeshTest, TransformIdentityAndRotation)
{
  std::vector<Vec3> points = { { 1., 0., 0. }, { 0., 1., 0. } };
  std::vector<real_t> ds = { 1., 1. };
  Vec3 cog = { 0., 0., 0. };

  TestRigidBody2DFromReference rb(points, ds, cog, 1., 1.);

  const auto& p0 = rb.mesh().world_points()[0];
  EXPECT_NEAR(p0[0], 1., 1e-12);
  EXPECT_NEAR(p0[1], 0., 1e-12);

  rb.set_state({ 0., 0., 0. }, { 0., 0., 0. }, std::acos(-1.) / 2., 0.);

  const auto& r0 = rb.mesh().world_points()[0];
  EXPECT_NEAR(r0[0], 0., 1e-12);
  EXPECT_NEAR(r0[1], 1., 1e-12);
  EXPECT_NEAR(r0[2], 0., 1e-12);
}

TEST(RigidBody2DIntegrationTest, UniformTractionSymmetry)
{
  std::vector<Vec3> points = { { 1., 0., 0. }, { -1., 0., 0. } };
  std::vector<real_t> ds = { 1., 1. };
  const Vec3 cog = { 0., 0., 0. };

  RigidBody2DFromReference rb(points, ds, cog, 1., 1.);

  std::vector<Vec3> traction = { { 2., 3., 0. }, { 2., 3., 0. } };

  const Vec3 F = rb.mesh().integrate_force(traction);
  const Vec3 T =
    rb.mesh().integrate_torque_about_com(traction, rb.center_of_mass());

  EXPECT_NEAR(F[0], 4., 1e-12);
  EXPECT_NEAR(F[1], 6., 1e-12);
  EXPECT_NEAR(F[2], 0., 1e-12);

  EXPECT_NEAR(T[0], 0., 1e-12);
  EXPECT_NEAR(T[1], 0., 1e-12);
  EXPECT_NEAR(T[2], 0., 1e-12);
}

TEST(RigidBody2DIntegrationTest, PureCoupleTorque)
{
  std::vector<Vec3> points = { { 1., 0., 0. }, { -1., 0., 0. } };
  std::vector<real_t> ds = { 1., 1. };
  const Vec3 cog = { 0., 0., 0. };

  RigidBody2DFromReference rb(points, ds, cog, 2., 4.);

  std::vector<Vec3> traction = { { 0., 1., 0. }, { 0., -1., 0. } };

  const Vec3 F = rb.mesh().integrate_force(traction);
  const Vec3 T =
    rb.mesh().integrate_torque_about_com(traction, rb.center_of_mass());

  EXPECT_NEAR(F[0], 0., 1e-12);
  EXPECT_NEAR(F[1], 0., 1e-12);
  EXPECT_NEAR(F[2], 0., 1e-12);

  EXPECT_NEAR(T[0], 0., 1e-12);
  EXPECT_NEAR(T[1], 0., 1e-12);
  EXPECT_NEAR(T[2], 2., 1e-12);
}

TEST(RigidBody2DDynamicsTest, ZeroTractionPreservesLinearAndAngularMomentum)
{
  std::vector<Vec3> points = { { 0., 0., 0. } };
  std::vector<real_t> ds = { 1. };

  TestRigidBody2DFromReference rb(points, ds, { 0., 0., 0. }, 2., 3.);
  rb.set_state({ 1., 2., 3. }, { 0.5, -0.25, 0.75 }, 0.2, 0.3);

  const Vec3 v0 = rb.com_velocity();
  const real_t Lz0 = rb.angular_momentum_z();

  rb.step(0.1, std::vector<Vec3>{ { 0., 0., 0. } });

  const Vec3 v1 = rb.com_velocity();
  const real_t Lz1 = rb.angular_momentum_z();

  EXPECT_NEAR(v1[0], v0[0], 1e-12);
  EXPECT_NEAR(v1[1], v0[1], 1e-12);
  EXPECT_NEAR(v1[2], v0[2], 1e-12);

  EXPECT_NEAR(Lz1, Lz0, 1e-12);
}

TEST(RigidBody2DDynamicsTest, ConstantTractionAcceleratesCOM)
{
  std::vector<Vec3> points = { { 0., 0., 0. } };
  std::vector<real_t> ds = { 2. };
  const Vec3 cog = { 0., 0., 0. };

  RigidBody2DFromReference rb(points, ds, cog, 4., 1.);

  rb.step(0.2, std::vector<Vec3>{ { 1., 0., 0. } });

  const Vec3 v = rb.com_velocity();
  const Vec3 x = rb.center_of_mass();

  EXPECT_NEAR(v[0], 0.1, 1e-12);
  EXPECT_NEAR(v[1], 0., 1e-12);
  EXPECT_NEAR(v[2], 0., 1e-12);

  EXPECT_NEAR(x[0], 0.02, 1e-12);
  EXPECT_NEAR(x[1], 0., 1e-12);
  EXPECT_NEAR(x[2], 0., 1e-12);
}

TEST(RigidBody2DContractTest, OrderingIsStable)
{
  std::vector<Vec3> points = { { 0., 0., 0. }, { 1., 0., 0. }, { 2., 0., 0. } };
  std::vector<real_t> ds = { 0.5, 1.0, 1.5 };

  RigidBody2DFromReference rb(points, ds, { 0., 0., 0. }, 1., 1.);

  const auto& ref = rb.mesh().reference_points();
  const auto& world = rb.mesh().world_points();

  ASSERT_EQ(ref.size(), world.size());
  for (size_t i = 0; i < ref.size(); ++i) {
    EXPECT_NEAR(ref[i][0], world[i][0], 1e-12);
    EXPECT_NEAR(ref[i][1], world[i][1], 1e-12);
    EXPECT_NEAR(ref[i][2], world[i][2], 1e-12);
  }
}

} // namespace
} // namespace Models
} // namespace ELFF
