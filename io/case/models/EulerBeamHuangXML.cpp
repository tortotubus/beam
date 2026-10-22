#include "elff/io/case/models/EulerBeamHuangXML.hpp"
#include "elff/io/case/InitialLineXML.hpp"
#include "elff/models/beam/EulerBeamHuangBuilder.hpp"

namespace ELFF::IO::Case {
using namespace Models;
namespace {
EulerBeam::EulerBeamBCType boundary_type(const CXX::XMLNode& boundary, const char* end) {
  const auto value = boundary.string(end);
  if (value == "free") return EulerBeam::free_bc;
  if (value == "simple" || value == "pinned") return EulerBeam::simple_bc;
  if (value == "clamped") return EulerBeam::clamped_bc;
  boundary.invalid("Unsupported EulerBeamHuang boundary condition.");
}
}
EulerBeamHuangDefinition EulerBeamHuangXML::parse(const CXX::XMLNode& model) {
  if (model.name() != tag) model.invalid("Expected <EulerBeamHuang>.");
  const auto geometry = model.child("Geometry");
  const auto material = model.child("Material");
  if (auto damping = optional_xml_child(material, "BendingDamping"))
    damping->invalid("BendingDamping is supported only by BendingRodSegments.");
  const bool v2 = bool(optional_xml_child(model, "InitialCondition"));
  EulerBeamHuangDefinition d;
  d.source_line = static_cast<long>(model.line());
  d.length = geometry.number("length");
  const int nodes = geometry.integer("nodes");
  if (nodes < 2) geometry.invalid("EulerBeamHuang requires at least two nodes.");
  d.nodes = static_cast<size_t>(nodes);
  d.bending_rigidity = v2 ? material.child("BendingRigidity").text_number() : material.number("bendingRigidity");
  d.linear_density = v2 ? material.child("LinearDensity").text_number() : material.number("linearDensity");
  d.boundary_conditions.end[0] = EulerBeam::left;
  d.boundary_conditions.end[1] = EulerBeam::right;
  if (v2) {
    const auto boundaries = model.children("EulerBeamBoundary");
    if (boundaries.size() != 2) model.invalid("Expected two EulerBeamBoundary elements.");
    bool seen[2] = {};
    for (const auto& boundary : boundaries) {
      const auto end = boundary.string("end");
      if (end != "left" && end != "right") boundary.invalid("Invalid beam boundary end.");
      const int i = end == "left" ? 0 : 1;
      if (seen[i]) boundary.invalid("Duplicate beam boundary end.");
      seen[i] = true;
      if (boundary.string("type") == "simple") boundary.invalid("Use pinned for a v2 simple support.");
      d.boundary_conditions.type[i] = boundary_type(boundary, "type");
    }
  } else {
    const auto boundary = model.child("BoundaryConditions");
    d.boundary_conditions.type[0] = boundary_type(boundary, "left");
    d.boundary_conditions.type[1] = boundary_type(boundary, "right");
  }
  if (const auto error = EulerBeamHuangBuilder::validation_error(d, InitialLineDefinition{}))
    model.invalid(error);
  return d;
}
InitialLineTemplate EulerBeamHuangXML::parse_placement(
    const CXX::XMLNode& model, const EulerBeamHuangDefinition& d,
    unsigned dimensions, bool allow_random) {
  auto placement = parse_initial_line_template(model, true, dimensions, allow_random);
  if (const auto error = EulerBeamHuangBuilder::validation_error(d, placement.fixed))
    initial_line_xml(model).invalid(error);
  return placement;
}
} // namespace ELFF::IO::Case
