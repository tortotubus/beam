#include "elff/c/models/beam/IBEulerBeamHuang.h"

#include "elff/c/models/beam/IBEulerBeamBCs.hpp"
#include "elff/config/config.hpp"
#include "elff/general/error.hpp"
#include "elff/models/beam/IBEulerBeamHuang.hpp"

#include <array>
#include <cmath>

using namespace ELFF::Models;
using namespace ELFF;

namespace {
EulerBeamMesh make_direction_initial_mesh(vertex_t s0,
                                          int nodes,
                                          double length,
                                          vertex_t direction);

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
  const vertex_t direction = { std::cos(theta), std::sin(theta), 0.0 };
  return make_direction_initial_mesh(s0, nodes, length, direction);
}

EulerBeamMesh
make_direction_initial_mesh(vertex_t s0,
                            int      nodes,
                            double   length,
                            vertex_t direction)
{
  const double direction_norm =
    std::sqrt(direction.x * direction.x + direction.y * direction.y +
              direction.z * direction.z);
  ELFF_VERIFY(std::isfinite(direction_norm) && direction_norm > 0.0,
              "make_direction_initial_mesh(): direction must be finite and "
              "nonzero.\n");

  const std::array<real_t, 3> tangent = {
    static_cast<real_t>(direction.x / direction_norm),
    static_cast<real_t>(direction.y / direction_norm),
    static_cast<real_t>(direction.z / direction_norm)
  };

  EulerBeamMesh ic_mesh(nodes, length);
  auto& ic_centerline = ic_mesh.get_centerline();
  auto& ic_slope = ic_mesh.get_slope();
  auto& ic_velocity = ic_mesh.get_centerline_velocity();
  auto& ic_s = ic_mesh.get_curvilinear_axis();

  for (size_t ni = 0; ni < static_cast<size_t>(nodes); ++ni) {
    const real_t s = ic_s[ni];

    ic_centerline[ni][0] = s0.x + (length - s) * tangent[0];
    ic_centerline[ni][1] = s0.y + (length - s) * tangent[1];
    ic_centerline[ni][2] = s0.z + (length - s) * tangent[2];

    ic_slope[ni][0] = -tangent[0];
    ic_slope[ni][1] = -tangent[1];
    ic_slope[ni][2] = -tangent[2];

    ic_velocity[ni][0] = 0.0;
    ic_velocity[ni][1] = 0.0;
    ic_velocity[ni][2] = 0.0;
  }

  return ic_mesh;
}
} // namespace

extern "C"
{

ib_euler_beam_huang_t
ib_euler_beam_huang_new(vertex_t s0,
                        ib_euler_beam_bcs_t bcs,
                        double   length,
                        double   EI,
                        double   mu,
                        int      nodes)
{
  EulerBeamMesh ic_mesh = make_offset_initial_mesh(s0, nodes, length);
  EulerBeam::EulerBeamBCs boundary_conditions = ELFF::C::to_cpp_beam_bcs(bcs);

  auto* beam = new IBEulerBeamHuang(static_cast<real_t>(length),
                                    static_cast<real_t>(EI),
                                    static_cast<real_t>(mu),
                                    static_cast<size_t>(nodes),
                                    boundary_conditions);

  beam->apply_initial_condition(ic_mesh);
  return beam;
}

ib_euler_beam_huang_t
ib_euler_beam_huang_new_theta(vertex_t s0,
                              ib_euler_beam_bcs_t bcs,
                              double   length,
                              double   EI,
                              double   mu,
                              int      nodes,
                              double   theta)
{
  EulerBeamMesh ic_mesh = make_theta_initial_mesh(s0, nodes, length, theta);
  EulerBeam::EulerBeamBCs boundary_conditions = ELFF::C::to_cpp_beam_bcs(bcs);

  auto* beam = new IBEulerBeamHuang(static_cast<real_t>(length),
                                    static_cast<real_t>(EI),
                                    static_cast<real_t>(mu),
                                    static_cast<size_t>(nodes),
                                    boundary_conditions);

  beam->apply_initial_condition(ic_mesh);
  return beam;
}

ib_euler_beam_huang_t
ib_euler_beam_huang_new_direction(vertex_t s0,
                                  ib_euler_beam_bcs_t bcs,
                                  double   length,
                                  double   EI,
                                  double   mu,
                                  int      nodes,
                                  vertex_t direction)
{
  EulerBeamMesh ic_mesh =
    make_direction_initial_mesh(s0, nodes, length, direction);
  EulerBeam::EulerBeamBCs boundary_conditions = ELFF::C::to_cpp_beam_bcs(bcs);

  auto* beam = new IBEulerBeamHuang(static_cast<real_t>(length),
                                    static_cast<real_t>(EI),
                                    static_cast<real_t>(mu),
                                    static_cast<size_t>(nodes),
                                    boundary_conditions);

  beam->apply_initial_condition(ic_mesh);
  return beam;
}

void
ib_euler_beam_huang_set_initial_velocity(ib_euler_beam_huang_t handle,
                                         const vertex_t*       velocity,
                                         int                   nodes)
{
  auto* beam = reinterpret_cast<IBEulerBeamHuang*>(handle);
  EulerBeamMesh initial_mesh = beam->get_mesh();
  auto& initial_velocity = initial_mesh.get_centerline_velocity();

  ELFF_VERIFY(velocity != nullptr,
              "ib_euler_beam_huang_set_initial_velocity(): velocity is null.\n");
  ELFF_VERIFY(nodes == static_cast<int>(initial_velocity.size()),
              "ib_euler_beam_huang_set_initial_velocity(): node count "
              "mismatch.\n");

  for (int i = 0; i < nodes; ++i) {
    initial_velocity[i][0] = velocity[i].x;
    initial_velocity[i][1] = velocity[i].y;
    initial_velocity[i][2] = velocity[i].z;
  }

  beam->apply_initial_condition(initial_mesh);
}

void
ib_euler_beam_huang_set_implicit_bending(ib_euler_beam_huang_t handle,
                                         int                   enabled)
{
  reinterpret_cast<IBEulerBeamHuang*>(handle)->set_implicit_bending(enabled != 0);
}

void
ib_euler_beam_huang_destroy(ib_euler_beam_huang_t handle)
{
  delete reinterpret_cast<IBEulerBeamHuang*>(handle);
}

} // extern "C"
