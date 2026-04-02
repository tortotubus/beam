#include "elff/c/models/beam/IBEulerBeamGGL.h"

#include "elff/config/config.hpp"
#include "elff/models/beam/IBEulerBeamGGL.hpp"

#include <cmath>

using namespace ELFF::Models;
using namespace ELFF;

namespace {
EulerBeam::EulerBeamBCs
make_beam_bcs(int bc_type_1, int bc_type_2)
{
  return {
    .end = { EulerBeam::left, EulerBeam::right },
    .type = { static_cast<EulerBeam::EulerBeamBCType>(bc_type_1),
              static_cast<EulerBeam::EulerBeamBCType>(bc_type_2) }
  };
}

EulerBeamMesh
make_offset_initial_mesh(vertex_t s0, int nodes, double length)
{
  EulerBeamMesh ic_mesh(nodes, length);
  auto& ic_centerline = ic_mesh.get_centerline();
  auto& ic_s = ic_mesh.get_curvilinear_axis();

  for (size_t ni = 0; ni < static_cast<size_t>(nodes); ++ni) {
    static_cast<void>(ic_s[ni]);
    ic_centerline[ni][0] += s0.x;
    ic_centerline[ni][1] += s0.y;
    ic_centerline[ni][2] += s0.z;
  }

  return ic_mesh;
}

EulerBeamMesh
make_theta_initial_mesh(vertex_t s0, int nodes, double length, double theta)
{
  const real_t x0 = s0.x;
  const real_t y0 = s0.y;
  const real_t z0 = s0.z;

  EulerBeamMesh ic_mesh(nodes, length);
  auto& ic_centerline = ic_mesh.get_centerline();
  auto& ic_slope = ic_mesh.get_slope();
  auto& ic_velocity = ic_mesh.get_centerline_velocity();
  auto& ic_s = ic_mesh.get_curvilinear_axis();

  for (size_t ni = 0; ni < static_cast<size_t>(nodes); ++ni) {
    const real_t s = ic_s[ni];

    ic_centerline[ni][0] = x0 + (length - s) * std::cos(theta);
    ic_centerline[ni][1] = y0 + (length - s) * std::sin(theta);
    ic_centerline[ni][2] = z0;

    ic_slope[ni][0] = -std::cos(theta);
    ic_slope[ni][1] = -std::sin(theta);
    ic_slope[ni][2] = 0.0;

    ic_velocity[ni][0] = 0.0;
    ic_velocity[ni][1] = 0.0;
    ic_velocity[ni][2] = 0.0;
  }

  return ic_mesh;
}
} // namespace

extern "C"
{

ib_euler_beam_ggl_t
ib_euler_beam_ggl_new(vertex_t s0,
                      int      bc_type_1,
                      int      bc_type_2,
                      double   length,
                      double   EI,
                      double   mu,
                      int      nodes,
                      double   r_penalty)
{
  EulerBeamMesh ic_mesh = make_offset_initial_mesh(s0, nodes, length);
  EulerBeam::EulerBeamBCs bcs = make_beam_bcs(bc_type_1, bc_type_2);

  auto* beam = new IBEulerBeamGGL(static_cast<real_t>(length),
                                  static_cast<real_t>(EI),
                                  static_cast<real_t>(mu),
                                  static_cast<size_t>(nodes),
                                  bcs,
                                  static_cast<real_t>(r_penalty));

  beam->apply_initial_condition(ic_mesh);
  return beam;
}

ib_euler_beam_ggl_t
ib_euler_beam_ggl_new_theta(vertex_t s0,
                            int      bc_type_1,
                            int      bc_type_2,
                            double   length,
                            double   EI,
                            double   mu,
                            int      nodes,
                            double   r_penalty,
                            double   theta)
{
  static_cast<void>(bc_type_1);
  static_cast<void>(bc_type_2);

  EulerBeamMesh ic_mesh = make_theta_initial_mesh(s0, nodes, length, theta);
  const real_t x0 = s0.x;
  const real_t y0 = s0.y;
  const real_t z0 = s0.z;

  EulerBeam::EulerBeamBCs boundary_conditions = {
    .end = { EulerBeam::left, EulerBeam::right },
    .type = { EulerBeam::free_bc, EulerBeam::simple_bc },
    .vals = { { .position = {} }, { .position = { x0, y0, z0 } } }
  };

  auto* beam = new IBEulerBeamGGL(static_cast<real_t>(length),
                                  static_cast<real_t>(EI),
                                  static_cast<real_t>(mu),
                                  static_cast<size_t>(nodes),
                                  boundary_conditions,
                                  static_cast<real_t>(r_penalty));

  beam->apply_initial_condition(ic_mesh);
  return beam;
}

void
ib_euler_beam_ggl_destroy(ib_euler_beam_ggl_t handle)
{
  delete reinterpret_cast<IBEulerBeamGGL*>(handle);
}

} // extern "C"
