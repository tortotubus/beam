#pragma once
#include <string>
namespace insertion_fixture {
inline std::string model(bool beam, const std::string& anchor, const std::string& direction,
                         const std::string& insertion = "") {
  const std::string tag = beam ? "EulerBeamHuang" : "BendingRodSegmentsImplicit";
  return "<" + tag + "><Geometry length='0.5' nodes='4'/>"
    "<Material bendingRigidity='0.01' linearDensity='1'/>"
    "<InitialLine anchorEnd='left'>" + anchor + direction + "</InitialLine>" +
    (beam ? "<BoundaryConditions left='free' right='free'/>" : "<RodBoundary type='free'/>") +
    insertion + "</" + tag + ">";
}
inline std::string fixed_anchor = "<Anchor x='1' y='1'/>";
inline std::string fixed_direction = "<Direction x='1' y='0'/>";
inline std::string random_anchor = "<RandomAnchor><Min x='1' y='1' z='0'/><Max x='9' y='9' z='0'/></RandomAnchor>";
inline std::string group = "<Insertion count='3' clearance='0.1' maxAttempts='100'/>";
inline std::string domain = "<InsertionDomain><Min x='0' y='0' z='0'/><Max x='10' y='10' z='0'/></InsertionDomain>";
inline std::string document(const std::string& models, const std::string& bounds = domain,
                            const char* policy = "round-robin", const char* dimensions = "2") {
  return std::string("<Case formatVersion='1'><ELFF schemaVersion='1' dimensions='") + dimensions +
    "' units='SI'>" + bounds + "<Models distribution='" + policy + "'>" + models + "</Models></ELFF></Case>";
}
inline std::string mixed(const char* policy = "round-robin") {
  return document(model(false, random_anchor, "<RandomDirection/>", group) +
                  model(true, "<Anchor x='5' y='5'/>", fixed_direction) +
                  model(true, random_anchor, "<RandomDirection/>", group), domain, policy);
}
inline std::string replace(std::string text, const std::string& from, const std::string& to) {
  text.replace(text.find(from), from.size(), to); return text;
}
}
