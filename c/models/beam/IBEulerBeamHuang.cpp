#include "elff/c/models/beam/IBEulerBeamHuang.h"

#include "elff/c/models/beam/IBEulerBeamBCs.hpp"
#include "elff/config/config.hpp"
#include "elff/models/beam/EulerBeamHuangBuilder.hpp"
#include "elff/models/ibm/IBRuntime.hpp"

#include <array>
#include <cmath>
#include <memory>

using namespace ELFF::Models;
using namespace ELFF;

extern "C"
{

ib_euler_beam_huang_t
ib_euler_beam_huang_new(ib_runtime_t runtime,
                        vertex_t s0,
                        ib_euler_beam_bcs_t bcs,
                        double   length,
                        double   EI,
                        double   mu,
                        int      nodes)
{
  if (!runtime)
    return nullptr;
  try {
  const std::array<real_t, 3> anchor = {s0.x, s0.y, s0.z};
  EulerBeamMesh ic_mesh = make_huang_offset_initial_mesh(
    anchor, static_cast<size_t>(nodes), static_cast<real_t>(length));
  EulerBeam::EulerBeamBCs boundary_conditions = ELFF::C::to_cpp_beam_bcs(bcs);

  auto beam = build_euler_beam_huang(std::move(ic_mesh), length, EI, mu,
                                     static_cast<size_t>(nodes),
                                     boundary_conditions);
  return reinterpret_cast<IBRuntime*>(runtime)->adopt_owned_model(
    IBRuntimeModelType::EulerBeamHuang, std::move(beam));
  } catch (...) {
    return nullptr;
  }
}

ib_euler_beam_huang_t
ib_euler_beam_huang_new_theta(ib_runtime_t runtime,
                              vertex_t s0,
                              ib_euler_beam_bcs_t bcs,
                              double   length,
                              double   EI,
                              double   mu,
                              int      nodes,
                              double   theta)
{
  const vertex_t direction = {std::cos(theta), std::sin(theta), 0.0};
  return ib_euler_beam_huang_new_direction(
    runtime, s0, bcs, length, EI, mu, nodes, direction);
}

ib_euler_beam_huang_t
ib_euler_beam_huang_new_direction(ib_runtime_t runtime,
                                  vertex_t s0,
                                  ib_euler_beam_bcs_t bcs,
                                  double   length,
                                  double   EI,
                                  double   mu,
                                  int      nodes,
                                  vertex_t direction)
{
  if (!runtime)
    return nullptr;
  try {
  const std::array<real_t, 3> anchor = {s0.x, s0.y, s0.z};
  const std::array<real_t, 3> tangent = {
    direction.x, direction.y, direction.z};
  EulerBeamMesh ic_mesh = make_huang_direction_initial_mesh(
    anchor, static_cast<size_t>(nodes), static_cast<real_t>(length), tangent);
  EulerBeam::EulerBeamBCs boundary_conditions = ELFF::C::to_cpp_beam_bcs(bcs);

  auto beam = build_euler_beam_huang(std::move(ic_mesh), length, EI, mu,
                                     static_cast<size_t>(nodes),
                                     boundary_conditions);
  return reinterpret_cast<IBRuntime*>(runtime)->adopt_owned_model(
    IBRuntimeModelType::EulerBeamHuang, std::move(beam));
  } catch (...) {
    return nullptr;
  }
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

} // extern "C"
