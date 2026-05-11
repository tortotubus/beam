#include <gtest/gtest.h>

#include "elff/c/models/beam/IBEulerBeamADDM.h"
#include "elff/c/models/beam/IBEulerBeamGGL.h"
#include "elff/c/models/beam/IBEulerBeamHuang.h"
#include "elff/c/models/beam/IBEulerBeamPenalty.h"
#include "elff/c/models/ibm/IBForceCoupled.h"

#define protected public
#include "elff/models/beam/IBEulerBeamADDM.hpp"
#undef protected

namespace {

ib_euler_beam_bcs_t
make_bcs(int left_type, int right_type)
{
  ib_euler_beam_bcs_t bcs = {};
  bcs.end[0] = IB_EULER_BEAM_BC_LEFT;
  bcs.end[1] = IB_EULER_BEAM_BC_RIGHT;
  bcs.type[0] = left_type;
  bcs.type[1] = right_type;
  return bcs;
}

void
expect_addm_constructs(ib_euler_beam_addm_t handle)
{
  ASSERT_NE(handle, nullptr);
  EXPECT_GT(ib_force_coupled_get_number_of_nodes(handle), 0);
  ib_euler_beam_addm_destroy(handle);
}

void
expect_penalty_constructs(ib_euler_beam_penalty_t handle)
{
  ASSERT_NE(handle, nullptr);
  EXPECT_GT(ib_force_coupled_get_number_of_nodes(handle), 0);
  ib_euler_beam_penalty_destroy(handle);
}

void
expect_ggl_constructs(ib_euler_beam_ggl_t handle)
{
  ASSERT_NE(handle, nullptr);
  EXPECT_GT(ib_force_coupled_get_number_of_nodes(handle), 0);
  ib_euler_beam_ggl_destroy(handle);
}

void
expect_huang_constructs(ib_euler_beam_huang_t handle)
{
  ASSERT_NE(handle, nullptr);
  EXPECT_GT(ib_force_coupled_get_number_of_nodes(handle), 0);
  ib_euler_beam_huang_destroy(handle);
}

} // namespace

TEST(IBEulerBeamCAPITest, ConstructsSolversWithFullBoundaryConditionStruct)
{
  const vertex_t s0 = { 0.0, 0.0, 0.0 };
  const double length = 1.0;
  const double EI = 1.0;
  const double mu = 1.0;
  const int nodes = 8;
  const double r_penalty = 1.0;

  ib_euler_beam_bcs_t free_simple =
    make_bcs(IB_EULER_BEAM_BC_FREE, IB_EULER_BEAM_BC_SIMPLE);

  ib_euler_beam_bcs_t simple_simple =
    make_bcs(IB_EULER_BEAM_BC_SIMPLE, IB_EULER_BEAM_BC_SIMPLE);
  simple_simple.vals[0].position = { 0.0, 0.0, 0.0 };
  simple_simple.vals[1].position = { 1.0, 0.0, 0.0 };

  ib_euler_beam_bcs_t clamped_free =
    make_bcs(IB_EULER_BEAM_BC_CLAMPED, IB_EULER_BEAM_BC_FREE);
  clamped_free.vals[0].position = { 0.0, 0.0, 0.0 };
  clamped_free.vals[0].slope = { 1.0, 0.0, 0.0 };

  expect_addm_constructs(ib_euler_beam_addm_new(
    s0, free_simple, length, EI, mu, nodes, r_penalty));
  expect_addm_constructs(ib_euler_beam_addm_new(
    s0, simple_simple, length, EI, mu, nodes, r_penalty));
  expect_addm_constructs(ib_euler_beam_addm_new(
    s0, clamped_free, length, EI, mu, nodes, r_penalty));

  expect_penalty_constructs(ib_euler_beam_penalty_new(
    s0, free_simple, length, EI, mu, nodes, r_penalty));
  expect_penalty_constructs(ib_euler_beam_penalty_new(
    s0, simple_simple, length, EI, mu, nodes, r_penalty));
  expect_penalty_constructs(ib_euler_beam_penalty_new(
    s0, clamped_free, length, EI, mu, nodes, r_penalty));

  expect_ggl_constructs(ib_euler_beam_ggl_new(
    s0, free_simple, length, EI, mu, nodes, r_penalty));
  expect_ggl_constructs(ib_euler_beam_ggl_new(
    s0, simple_simple, length, EI, mu, nodes, r_penalty));
  expect_ggl_constructs(ib_euler_beam_ggl_new(
    s0, clamped_free, length, EI, mu, nodes, r_penalty));

  expect_huang_constructs(
    ib_euler_beam_huang_new(s0, free_simple, length, EI, mu, nodes));
  expect_huang_constructs(
    ib_euler_beam_huang_new(s0, simple_simple, length, EI, mu, nodes));
  expect_huang_constructs(
    ib_euler_beam_huang_new(s0, clamped_free, length, EI, mu, nodes));
}

TEST(IBEulerBeamCAPITest, ThetaConstructorUsesSuppliedBoundaryConditions)
{
  const vertex_t s0 = { 2.0, 3.0, 4.0 };
  ib_euler_beam_bcs_t bcs =
    make_bcs(IB_EULER_BEAM_BC_CLAMPED, IB_EULER_BEAM_BC_FREE);
  bcs.vals[0].position = { 2.0, 3.0, 4.0 };
  bcs.vals[0].slope = { 0.5, 0.25, 0.0 };
  bcs.vals[1].force = { 7.0, 8.0, 9.0 };

  ib_euler_beam_addm_t handle =
    ib_euler_beam_addm_new_theta(s0, bcs, 1.0, 1.0, 1.0, 8, 1.0, 0.3);
  ASSERT_NE(handle, nullptr);

  auto* beam = reinterpret_cast<ELFF::Models::IBEulerBeamADDM*>(handle);
  EXPECT_EQ(beam->boundary_conditions.end[0], ELFF::Models::EulerBeam::left);
  EXPECT_EQ(beam->boundary_conditions.end[1], ELFF::Models::EulerBeam::right);
  EXPECT_EQ(beam->boundary_conditions.type[0],
            ELFF::Models::EulerBeam::clamped_bc);
  EXPECT_EQ(beam->boundary_conditions.type[1],
            ELFF::Models::EulerBeam::free_bc);
  EXPECT_DOUBLE_EQ(beam->boundary_conditions.vals[0].position[1], 3.0);
  EXPECT_DOUBLE_EQ(beam->boundary_conditions.vals[0].slope[0], 0.5);
  EXPECT_DOUBLE_EQ(beam->boundary_conditions.vals[1].force[2], 9.0);

  ib_euler_beam_addm_destroy(handle);
}
