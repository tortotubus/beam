#include "elff/models/beam/EulerBeamHuangBuilder.hpp"
#include <climits>

namespace ELFF::Models {
const char* EulerBeamHuangBuilder::validation_error(
    const EulerBeamHuangDefinition& d, const ELFF::IO::Case::InitialLineDefinition& p) noexcept {
  if (!std::isfinite(d.length) || d.length <= 0 ||
      !std::isfinite(d.bending_rigidity) || d.bending_rigidity <= 0 ||
      !std::isfinite(d.linear_density) || d.linear_density <= 0 ||
      d.nodes < 2 || d.nodes > INT_MAX)
    return "EulerBeamHuang length, material values, and node count must be positive (nodes >= 2 and <= INT_MAX).";
  const auto norm = std::hypot(p.direction[0], p.direction[1], p.direction[2]);
  if (!std::isfinite(norm) || norm == 0)
    return "EulerBeamHuang direction must be finite and nonzero.";
  for (size_t k = 0; k < 3; ++k)
    if (!std::isfinite(p.anchor[k]) ||
        !std::isfinite(p.anchor[k] + d.length * (p.direction[k] / norm)))
      return "EulerBeamHuang initial positions are outside the configured real range.";
  for (auto type : d.boundary_conditions.type)
    if (type != EulerBeam::free_bc && type != EulerBeam::simple_bc &&
        type != EulerBeam::clamped_bc)
      return "Unsupported EulerBeamHuang boundary condition.";
  return nullptr;
}

std::unique_ptr<IBEulerBeamHuang> EulerBeamHuangBuilder::build(
    const EulerBeamHuangDefinition& d, const ELFF::IO::Case::InitialLineDefinition& p) {
  const auto error = validation_error(d, p);
  ELFF_VERIFY(!error, (error ? error : "Invalid Huang definition"));
  auto mesh = make_huang_direction_initial_mesh(p.anchor, d.nodes, d.length,
                                               p.direction, p.anchor_at_right);
  auto bcs = d.boundary_conditions;
  bcs.end[0] = EulerBeam::left;
  bcs.end[1] = EulerBeam::right;
  for (size_t end = 0; end < 2; ++end) {
    const auto node = end == 0 ? 0 : d.nodes - 1;
    for (size_t k = 0; k < 3; ++k) {
      bcs.vals[end].position[k] = mesh.get_centerline(node)[k];
      bcs.vals[end].slope[k] = mesh.get_slope(node)[k];
    }
  }
  return build_euler_beam_huang(std::move(mesh), d.length, d.bending_rigidity,
                               d.linear_density, d.nodes, bcs);
}
} // namespace ELFF::Models
