#include "elff/models/rigidbody/RigidBody2DModel.hpp"

#include <gtest/gtest.h>

#include <vector>

namespace ELFF {
namespace Models {
namespace {

using Vec3 = RigidBodyModel::Vec3;

class TestRigidBody2D : public RigidBody2DModel
{
public:
  TestRigidBody2D(const std::vector<Vec3>& points_ref,
                  const std::vector<real_t>& ds,
                  const Vec3& cog_ref,
                  real_t mass,
                  real_t Izz)
    : points_ref_(points_ref)
    , ds_(ds)
    , cog_ref_(cog_ref)
    , mass_(mass)
    , Izz_(Izz)
  {
    initialize();
  }

  void set_state(const Vec3& x_com_world,
                 const Vec3& v_com_world,
                 real_t theta,
                 real_t omega_z)
  {
    set_initial_state(x_com_world, v_com_world, theta, omega_z);
  }

protected:
  void define_reference_configuration(
    std::vector<Vec3>& points_ref,
    std::vector<real_t>& ds,
    Vec3& cog_ref,
    std::vector<Vec3>& normals_ref) const override
  {
    points_ref = points_ref_;
    ds = ds_;
    cog_ref = cog_ref_;
    normals_ref.clear();
  }

  void define_mass_properties(real_t& mass, real_t& Izz) const override
  {
    mass = mass_;
    Izz = Izz_;
  }

private:
  std::vector<Vec3> points_ref_;
  std::vector<real_t> ds_;
  Vec3 cog_ref_;
  real_t mass_ = 1.;
  real_t Izz_ = 1.;
};

TEST(RigidBody2DModelTest, ZeroTractionPreservesPlanarMomenta)
{
  TestRigidBody2D rb(std::vector<Vec3>{ { 0., 0., 0. } },
                     std::vector<real_t>{ 1. },
                     { 0., 0., 0. },
                     2.,
                     3.);
  rb.set_state({ 1., 2., 7. }, { 0.5, -0.25, 8. }, 0.2, 0.3);

  const Vec3 v0 = rb.com_velocity();
  const real_t Lz0 = rb.angular_momentum_z();

  rb.step(0.1, std::vector<Vec3>{ { 0., 0., 0. } });

  const Vec3 v1 = rb.com_velocity();
  EXPECT_NEAR(v1[0], v0[0], 1e-12);
  EXPECT_NEAR(v1[1], v0[1], 1e-12);
  EXPECT_NEAR(v1[2], 0., 1e-12);
  EXPECT_NEAR(rb.angular_momentum_z(), Lz0, 1e-12);
}

TEST(RigidBody2DModelTest, ConstantTractionAcceleratesCOM)
{
  TestRigidBody2D rb(std::vector<Vec3>{ { 0., 0., 0. } },
                     std::vector<real_t>{ 2. },
                     { 0., 0., 0. },
                     4.,
                     1.);

  rb.step(0.2, std::vector<Vec3>{ { 1., 0., 7. } });

  const Vec3 v = rb.com_velocity();
  const Vec3 x = rb.center_of_mass();
  EXPECT_NEAR(v[0], 0.1, 1e-12);
  EXPECT_NEAR(v[1], 0., 1e-12);
  EXPECT_NEAR(v[2], 0., 1e-12);

  EXPECT_NEAR(x[0], 0.02, 1e-12);
  EXPECT_NEAR(x[1], 0., 1e-12);
  EXPECT_NEAR(x[2], 0., 1e-12);
}

TEST(RigidBody2DModelTest, PureCoupleUpdatesAngularState)
{
  TestRigidBody2D rb(std::vector<Vec3>{ { 1., 0., 0. }, { -1., 0., 0. } },
                     std::vector<real_t>{ 1., 1. },
                     { 0., 0., 0. },
                     1.,
                     2.);
  rb.set_state({ 0., 0., 0. }, { 0., 0., 0. }, 0., 0.);

  rb.step(0.1, std::vector<Vec3>{ { 0., 1., 0. }, { 0., -1., 0. } });

  EXPECT_NEAR(rb.angular_velocity_z(), 0.1, 1e-12);
  EXPECT_NEAR(rb.angular_momentum_z(), 0.2, 1e-12);
  EXPECT_NEAR(rb.angle(), 0.01, 1e-12);
}

} // namespace
} // namespace Models
} // namespace ELFF
