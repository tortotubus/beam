#pragma once

#include "elff/models/beam/EulerBeamHuangDefinition.hpp"
#include "elff/io/case/ModelPlacement.hpp"
#include "elff/config/config.hpp"
#include "elff/general/error.hpp"
#include "elff/models/beam/IBEulerBeamHuang.hpp"

#include <array>
#include <cmath>
#include <memory>

namespace ELFF {
namespace Models {

/** @brief Construct initialized Huang models from owned settings and resolved placement. */
class EulerBeamHuangBuilder {
public:
  /** @brief Check settings and derived endpoints without constructing a model.
   * @return A diagnostic on failure, or nullptr for valid inputs. */
  static const char* validation_error(const EulerBeamHuangDefinition&,
                                     const ELFF::IO::Case::InitialLineDefinition&) noexcept;
  /** @brief Build the straight initial mesh and set boundary positions and slopes. */
  static std::unique_ptr<IBEulerBeamHuang> build(
      const EulerBeamHuangDefinition&, const ELFF::IO::Case::InitialLineDefinition&);
};

inline EulerBeamMesh
make_huang_offset_initial_mesh(const std::array<real_t, 3>& anchor,
                               size_t nodes,
                               real_t length)
{
  EulerBeamMesh mesh(nodes, length);
  auto& centerline = mesh.get_centerline();
  for (auto& point : centerline) {
    point[0] += anchor[0];
    point[1] += anchor[1];
    point[2] += anchor[2];
  }
  return mesh;
}

inline EulerBeamMesh
make_huang_direction_initial_mesh(const std::array<real_t, 3>& anchor,
                                  size_t nodes,
                                  real_t length,
                                  const std::array<real_t, 3>& direction,
                                  bool anchor_at_right = true)
{
  const real_t norm = std::hypot(direction[0], direction[1], direction[2]);
  ELFF_VERIFY(std::isfinite(norm) && norm > 0.0,
              "EulerBeamHuang direction must be finite and nonzero.\n");
  const std::array<real_t, 3> tangent = {
    direction[0] / norm, direction[1] / norm, direction[2] / norm};

  EulerBeamMesh mesh(nodes, length);
  auto& centerline = mesh.get_centerline();
  auto& slope = mesh.get_slope();
  auto& velocity = mesh.get_centerline_velocity();
  auto& s_axis = mesh.get_curvilinear_axis();
  for (size_t i = 0; i < nodes; ++i) {
    const real_t distance = anchor_at_right ? length - s_axis[i] : s_axis[i];
    const real_t slope_sign = anchor_at_right ? -1.0 : 1.0;
    for (size_t d = 0; d < 3; ++d) {
      centerline[i][d] = anchor[d] + distance * tangent[d];
      slope[i][d] = slope_sign * tangent[d];
      velocity[i][d] = 0.0;
    }
  }
  return mesh;
}

inline std::unique_ptr<IBEulerBeamHuang>
build_euler_beam_huang(EulerBeamMesh initial_mesh,
                       real_t length,
                       real_t EI,
                       real_t mu,
                       size_t nodes,
                       const EulerBeam::EulerBeamBCs& boundary_conditions)
{
  auto beam = std::make_unique<IBEulerBeamHuang>(
    length, EI, mu, nodes, boundary_conditions);
  beam->apply_initial_condition(initial_mesh);
  return beam;
}

} // namespace Models
} // namespace ELFF
