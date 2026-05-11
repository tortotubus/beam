#include "elff/c/models/beam/IBEulerBeamADDM.h"

#include "elff/c/models/beam/IBEulerBeamBCs.hpp"
#include "elff/config/config.hpp"
#include "elff/models/beam/IBEulerBeamADDM.hpp"

#include <cmath>

using namespace ELFF::Models;
using namespace ELFF;

namespace {
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

ib_euler_beam_addm_t
ib_euler_beam_addm_new(vertex_t s0,
                       ib_euler_beam_bcs_t bcs,
                       double   length,
                       double   EI,
                       double   mu,
                       int      nodes,
                       double   r_penalty)
{
  EulerBeamMesh ic_mesh = make_offset_initial_mesh(s0, nodes, length);
  EulerBeam::EulerBeamBCs boundary_conditions = ELFF::C::to_cpp_beam_bcs(bcs);

  auto* beam = new IBEulerBeamADDM(static_cast<real_t>(length),
                                   static_cast<real_t>(EI),
                                   static_cast<real_t>(mu),
                                   static_cast<size_t>(nodes),
                                   boundary_conditions,
                                   static_cast<real_t>(r_penalty));

  beam->apply_initial_condition(ic_mesh);
  return beam;
}

ib_euler_beam_addm_t
ib_euler_beam_addm_new_theta(vertex_t s0,
                             ib_euler_beam_bcs_t bcs,
                             double   length,
                             double   EI,
                             double   mu,
                             int      nodes,
                             double   r_penalty,
                             double   theta)
{
  EulerBeamMesh ic_mesh = make_theta_initial_mesh(s0, nodes, length, theta);
  EulerBeam::EulerBeamBCs boundary_conditions = ELFF::C::to_cpp_beam_bcs(bcs);

  auto* beam = new IBEulerBeamADDM(static_cast<real_t>(length),
                                   static_cast<real_t>(EI),
                                   static_cast<real_t>(mu),
                                   static_cast<size_t>(nodes),
                                   boundary_conditions,
                                   static_cast<real_t>(r_penalty));

  beam->apply_initial_condition(ic_mesh);
  return beam;
}

void
ib_euler_beam_addm_destroy(ib_euler_beam_addm_t handle)
{
  delete reinterpret_cast<IBEulerBeamADDM*>(handle);
}

} // extern "C"
