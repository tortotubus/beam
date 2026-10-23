#include "elff/io/case/models/BendingRodImplicitXML.hpp"
#include "elff/io/CXX/xmlCase.hpp"
#include "elff/io/case/InitialLineXML.hpp"
#include <climits>
#include <cmath>
#include <optional>

namespace ELFF::IO::Case {
using namespace Models;
namespace {
using IO::CXX::XMLNode;
using IO::CXX::XMLError;
double number_or(const XMLNode &node, const char *name, double fallback) {
  try { return node.number(name); }
  catch (const XMLError &error) {
    if (error.status() != ELFF_XML_NOT_FOUND) throw;
    return fallback;
  }
}
bool representable(double value) {
  return std::isfinite(value) && std::isfinite(static_cast<real_t>(value));
}
}

BendingRodImplicitDefinition parse_bending_rod_implicit_xml(elff_xml_node_t raw)
{
  return parse_bending_rod_implicit_xml(raw, IO::Case::parse_fixed_initial_line(IO::CXX::XMLNode(raw), false));
}

BendingRodImplicitDefinition parse_bending_rod_implicit_xml(elff_xml_node_t raw, const InitialLineDefinition& placement)
{
  const XMLNode model(raw);
  const bool v2 = model.name() == "BendingRodSegments";
  if (!v2 && model.name() != "BendingRodSegmentsImplicit")
    model.invalid("Expected <BendingRodSegments>.");
  BendingRodImplicitDefinition result;
  result.anchor_at_right = placement.anchor_at_right;
  result.source_line = static_cast<long>(model.line());
  auto &config = result.config;
  const auto geometry = model.child("Geometry");
  const int nodes = geometry.integer("nodes");
  config.length = geometry.number("length");
  if (nodes < 2 || nodes > INT_MAX / 6 || !representable(config.length) ||
      !(static_cast<real_t>(config.length) > 0))
    geometry.invalid("Rod requires positive finite length and 2 <= nodes <= INT_MAX/6.");
  config.links = nodes - 1;
  const auto material = model.child("Material");
  config.ei = v2 ? material.child("BendingRigidity").text_number() : material.number("bendingRigidity");
  const double density = v2 ? material.child("LinearDensity").text_number() : material.number("linearDensity");
  config.total_mass = config.length * density;
  const double h = config.length / config.links;
  if (auto damping = optional_xml_child(material, "BendingDamping")) {
    config.damping_rigidity = damping->text_number();
    if (!representable(config.damping_rigidity) || config.damping_rigidity < 0 ||
        !representable(config.damping_rigidity / h))
      damping->invalid("BendingDamping and derived joint damping must be finite and nonnegative.");
  }
  const double inertia = (config.total_mass / config.links) * h * h / 12;
  if (!representable(config.ei) || config.ei < 0 || !representable(density) || density <= 0 ||
      !representable(config.total_mass) || !(static_cast<real_t>(config.total_mass) > 0) ||
      !representable(h) || !(static_cast<real_t>(h) > 0) ||
      !representable(inertia) || !(static_cast<real_t>(inertia) > 0) ||
      !representable(config.ei / h))
    material.invalid("Rod material and derived mass/inertia must be finite and positive; bending rigidity may be zero.");
  const auto line = initial_line_xml(model);
  config.pin = Vec3(placement.anchor[0], placement.anchor[1], placement.anchor[2]);
  result.direction = Vec3(placement.direction[0], placement.direction[1], placement.direction[2]);
  const double norm = result.direction.stableNorm();
  if (!std::isfinite(norm) || norm <= 0)
    line.invalid("Rod direction must be finite and nonzero.");
  result.direction /= norm;
  const Vec3 end = config.pin + config.length * result.direction;
  for (int i = 0; i < 3; ++i)
    if (!representable(config.pin[i]) || !representable(end[i]))
      line.invalid("Rod initial positions are outside the configured real range.");
  if (v2) {
    const auto boundaries = model.children("BendingRodBoundary");
    if (boundaries.size() != 2) model.invalid("Expected two BendingRodBoundary elements.");
    bool seen[2] = {};
    int pins = 0;
    config.boundary = Boundary::Free;
    for (const auto& boundary : boundaries) {
      const auto end = boundary.string("end");
      if (end != "left" && end != "right") boundary.invalid("Invalid rod boundary end.");
      const int i = end == "left" ? 0 : 1;
      if (seen[i]) boundary.invalid("Duplicate rod boundary end.");
      seen[i] = true;
      const auto type = boundary.string("type");
      if (type != "free" && type != "pinned")
        boundary.invalid("Rod boundary must be pinned or free.");
      if (type == "pinned") {
        if (++pins > 1) boundary.invalid("Rod supports at most one pinned endpoint.");
        config.boundary = i == 0 ? Boundary::Pinned : Boundary::PinnedRight;
      }
    }
  } else {
    const auto boundary = model.child("RodBoundary");
    const auto type = boundary.string("type");
    if (type == "pinned") config.boundary = Boundary::Pinned;
    else if (type == "free") config.boundary = Boundary::Free;
    else boundary.invalid("Rod boundary must be 'pinned' or 'free'.");
  }
  if (auto gravity = optional_xml_child(model, "Gravity"))
    gravity->invalid("Gravity is no longer supported; supply external nodal forces.");
  if (auto solver = optional_xml_child(model, "ImplicitSolver")) {
    config.spectral_radius = number_or(*solver, "spectralRadius", config.spectral_radius);
    config.nonlinear_tolerance = number_or(*solver, "nonlinearTolerance", config.nonlinear_tolerance);
    try { config.max_iterations = solver->integer("maxIterations"); }
    catch (const XMLError &error) { if (error.status() != ELFF_XML_NOT_FOUND) throw; }
    if (!representable(config.spectral_radius) || config.spectral_radius < 0 || config.spectral_radius > 1 ||
        !representable(config.nonlinear_tolerance) || !(static_cast<real_t>(config.nonlinear_tolerance) > 0) ||
        config.max_iterations < 1)
      solver->invalid("Invalid implicit solver controls: radius must be in [0,1], tolerance and iterations positive.");
  }
  return result;
}
} // namespace ELFF::IO::Case
