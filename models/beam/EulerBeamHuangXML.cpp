#include "elff/models/beam/EulerBeamHuangXML.hpp"
#include "elff/io/CXX/xmlCase.hpp"
#include <cmath>
#include <string>

namespace ELFF::Models {
namespace {
bool parse_bc_type(const std::string& name, EulerBeam::EulerBeamBCType& type)
{
  if (name == "free")
    type = EulerBeam::free_bc;
  else if (name == "simple")
    type = EulerBeam::simple_bc;
  else if (name == "clamped")
    type = EulerBeam::clamped_bc;
  else
    return false;
  return true;
}

} // namespace

EulerBeamHuangDefinition parse_euler_beam_huang_xml(elff_xml_node_t node)
{
  const IO::CXX::XMLNode model(node);
  if (model.name() != "EulerBeamHuang")
    model.invalid("Expected <EulerBeamHuang>.");
  const auto geometry = model.child("Geometry");
  const auto material = model.child("Material");
  const auto initial_line = model.child("InitialLine");
  const auto anchor = initial_line.child("Anchor");
  const auto direction = initial_line.child("Direction");
  const auto boundary = model.child("BoundaryConditions");

  const double length = geometry.number("length");
  const int nodes = geometry.integer("nodes");
  const double bending = material.number("bendingRigidity");
  const double density = material.number("linearDensity");
  const auto anchor_value = anchor.xyz();
  const auto direction_value = direction.xyz();
  const auto anchor_end = initial_line.string("anchorEnd");
  const auto left_bc = boundary.string("left");
  const auto right_bc = boundary.string("right");
  const size_t line = model.line();
  EulerBeamHuangDefinition definition{};
  const bool finite_scalars =
    std::isfinite(length) && std::isfinite(bending) &&
    std::isfinite(density);
  bool finite_vectors = true;
  for (size_t d = 0; d < 3; ++d) {
    finite_vectors = finite_vectors && std::isfinite(anchor_value[d]) &&
                     std::isfinite(direction_value[d]);
  }
  if (!finite_scalars || !(length > 0.0) || !(bending > 0.0) ||
      !(density > 0.0) || nodes < 2) {
    model.invalid("EulerBeamHuang length, material values, and node count must "
              "be positive (nodes >= 2).");
  }
  if (!finite_vectors) {
    model.invalid("EulerBeamHuang anchor and direction values must be finite.");
  }
  const double direction_norm =
    std::sqrt(direction_value[0] * direction_value[0] +
              direction_value[1] * direction_value[1] +
              direction_value[2] * direction_value[2]);
  if (!std::isfinite(direction_norm) || direction_norm == 0.0) {
    model.invalid("EulerBeamHuang direction must be finite and nonzero.");
  }
  if (anchor_end != "left" && anchor_end != "right") {
    model.invalid("EulerBeamHuang anchorEnd must be 'left' or 'right'.");
  }

  EulerBeam::EulerBeamBCs bcs = {};
  bcs.end[0] = EulerBeam::left;
  bcs.end[1] = EulerBeam::right;
  if (!parse_bc_type(left_bc, bcs.type[0]) ||
      !parse_bc_type(right_bc, bcs.type[1])) {
    model.invalid("Unsupported EulerBeamHuang boundary condition.");
  }

  const bool anchor_at_right = anchor_end == "right";
  const real_t sign = anchor_at_right ? -1.0 : 1.0;
  for (size_t d = 0; d < 3; ++d) {
    const real_t a = static_cast<real_t>(anchor_value[d]);
    const real_t tangent = static_cast<real_t>(direction_value[d] /
                                               direction_norm);
    if (!std::isfinite(a) || !std::isfinite(tangent)) {
      model.invalid("EulerBeamHuang anchor or direction is outside the "
                "configured real range.");
    }
    const real_t left_position =
      anchor_at_right ? a + static_cast<real_t>(length) * tangent : a;
    const real_t right_position =
      anchor_at_right ? a : a + static_cast<real_t>(length) * tangent;
    bcs.vals[0].position[d] = left_position;
    bcs.vals[1].position[d] = right_position;
    bcs.vals[0].slope[d] = sign * tangent;
    bcs.vals[1].slope[d] = sign * tangent;
    definition.anchor[d] = a;
    definition.direction[d] = static_cast<real_t>(direction_value[d]);
  }

  definition.source_line = static_cast<long>(line);
  definition.length = static_cast<real_t>(length);
  definition.bending_rigidity = static_cast<real_t>(bending);
  definition.linear_density = static_cast<real_t>(density);
  if (!std::isfinite(definition.length) ||
      !std::isfinite(definition.bending_rigidity) ||
      !std::isfinite(definition.linear_density)) {
    model.invalid("EulerBeamHuang values are outside the configured real range.");
  }
  definition.nodes = static_cast<size_t>(nodes);
  definition.anchor_at_right = anchor_at_right;
  definition.boundary_conditions = bcs;
  return definition;
}

} // namespace ELFF::Models
