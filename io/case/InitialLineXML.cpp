#include "elff/io/case/InitialLineXML.hpp"
#include <algorithm>
#include <cmath>
#include <limits>

namespace ELFF::IO::Case {
using IO::CXX::XMLNode;
using IO::CXX::XMLError;
std::optional<XMLNode> optional_xml_child(const XMLNode& node, const char* name) {
  try { return node.child(name); }
  catch (const XMLError& error) {
    if (error.status() != ELFF_XML_NOT_FOUND) throw;
    return std::nullopt;
  }
}
XMLNode initial_line_xml(const XMLNode& model) {
  if (auto initial = optional_xml_child(model, "InitialCondition"))
    return initial->child("Line");
  return model.child("InitialLine");
}
namespace {
std::array<real_t, 3> vector(const XMLNode& node) {
  const auto values = node.xyz();
  std::array<real_t, 3> result{};
  for (size_t d = 0; d < 3; ++d) {
    result[d] = values[d];
    if (!std::isfinite(result[d])) node.invalid("Placement coordinates must be finite and representable.");
  }
  return result;
}
}
PlacementBox parse_placement_box(const XMLNode& node, unsigned dimensions, bool domain) {
  PlacementBox box;
  box.source_line = node.line();
  auto lo = node.child("Min"), hi = node.child("Max");
  if (dimensions == 3) { lo.number("z"); hi.number("z"); }
  box.min = vector(lo); box.max = vector(hi);
  for (unsigned d = 0; d < 3; ++d) {
    const real_t extent = box.max[d] - box.min[d];
    if (!std::isfinite(extent) || extent < 0 || (domain && d < dimensions && extent == 0))
      node.invalid("Placement bounds require finite ordered extents; domain extents must be positive.");
  }
  if (dimensions == 2 && (box.min[2] != 0 || box.max[2] != 0))
    node.invalid("Two-dimensional placement must lie in the XY plane.");
  return box;
}
InitialLineTemplate parse_initial_line_template(const XMLNode& model, bool huang,
                                                unsigned dimensions, bool allow_random) {
  auto line = initial_line_xml(model);
  InitialLineTemplate result;
  result.source_line = line.line();
  std::string end = "left";
  try { end = line.string("anchorEnd"); }
  catch (const XMLError& error) {
    if (huang || error.status() != ELFF_XML_NOT_FOUND) throw;
  }
  if (end != "left" && !(huang && end == "right"))
    line.invalid("Invalid anchorEnd: Huang requires left or right; rods allow only left or omission.");
  result.fixed.anchor_at_right = end == "right";
  auto anchor = optional_xml_child(line, "Anchor");
  auto random_anchor = optional_xml_child(line, "RandomAnchor");
  auto direction = optional_xml_child(line, "Direction");
  auto random_direction = optional_xml_child(line, "RandomDirection");
  if (bool(anchor) == bool(random_anchor) || bool(direction) == bool(random_direction))
    line.invalid("InitialLine requires exactly one fixed or random anchor and direction.");
  if (!allow_random && (random_anchor || random_direction))
    line.invalid("Random placement requires an insertion template.");
  if (anchor) result.fixed.anchor = vector(*anchor);
  else {
    result.anchors = parse_placement_box(*random_anchor, dimensions, false);
    result.fixed.anchor = result.anchors->min;
  }
  if (direction) {
    result.fixed.direction = vector(*direction);
    auto& v = result.fixed.direction;
    const real_t norm = std::hypot(v[0], v[1], v[2]);
    if (!std::isfinite(norm) || norm == 0) direction->invalid("InitialLine direction must be finite and nonzero.");
    for (auto& x : v) x /= norm;
  } else result.random_direction = true;
  if (dimensions == 2 && (result.fixed.anchor[2] != 0 || result.fixed.direction[2] != 0))
    line.invalid("Two-dimensional placement must lie in the XY plane.");
  return result;
}
InitialLineDefinition parse_fixed_initial_line(const XMLNode& model, bool huang) {
  // Existing non-insertion cases retain their unrestricted 3D placement.
  return parse_initial_line_template(model, huang, 3, false).fixed;
}
} // namespace ELFF::IO::Case
