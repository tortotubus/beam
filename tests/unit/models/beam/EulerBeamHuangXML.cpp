#include <gtest/gtest.h>
#include "elff/config/_config.hpp"
#include "elff/models/beam/EulerBeamHuangXML.hpp"
#include "elff/io/CXX/xmlCase.hpp"
#include <memory>
#include <string>

using namespace ELFF::Models;
using ELFF::IO::CXX::XMLError;

#ifdef ELFF_USE_LIBXML2
namespace {
const std::string valid = R"(<EulerBeamHuang>
<Geometry length="2" nodes="5"/>
<Material bendingRigidity="0.01" linearDensity="1.5"/>
<InitialLine anchorEnd="right"><Anchor x="1" y="2"/><Direction x="0" y="2"/></InitialLine>
<BoundaryConditions left="free" right="clamped"/>
</EulerBeamHuang>)";

EulerBeamHuangDefinition parse(const std::string& xml)
{
  elff_xml_document_t raw = nullptr;
  elff_xml_error error{};
  if (elff_xml_document_load_memory(xml.data(), xml.size(), &raw, &error) != ELFF_XML_OK)
    throw std::runtime_error(error.message);
  std::unique_ptr<elff_xml_document_handle, decltype(&elff_xml_document_destroy)>
    document(raw, elff_xml_document_destroy);
  elff_xml_node_t root = nullptr;
  if (elff_xml_document_root(raw, "EulerBeamHuang", &root, &error) != ELFF_XML_OK)
    throw std::runtime_error(error.message);
  return parse_euler_beam_huang_xml(root);
}
std::string changed(const std::string& from, const std::string& to)
{
  auto xml = valid;
  xml.replace(xml.find(from), from.size(), to);
  return xml;
}
}

TEST(EulerBeamHuangXMLTest, ParsesWithoutRuntimeAndSurvivesDocumentDestruction)
{
  const auto definition = parse(valid);
  EXPECT_EQ(definition.nodes, 5u);
  EXPECT_EQ(definition.source_line, 1);
  EXPECT_DOUBLE_EQ(definition.length, 2);
  EXPECT_TRUE(definition.anchor_at_right);
  EXPECT_DOUBLE_EQ(definition.boundary_conditions.vals[0].position[1], 4);
  EXPECT_DOUBLE_EQ(definition.boundary_conditions.vals[1].position[1], 2);
  EXPECT_DOUBLE_EQ(definition.boundary_conditions.vals[0].slope[1], -1);
  EXPECT_EQ(definition.boundary_conditions.type[1], EulerBeam::clamped_bc);
}

TEST(EulerBeamHuangXMLTest, SupportsLeftAnchoring)
{
  const auto definition = parse(changed("anchorEnd=\"right\"", "anchorEnd=\"left\""));
  EXPECT_FALSE(definition.anchor_at_right);
  EXPECT_DOUBLE_EQ(definition.boundary_conditions.vals[0].position[1], 2);
  EXPECT_DOUBLE_EQ(definition.boundary_conditions.vals[1].position[1], 4);
  EXPECT_DOUBLE_EQ(definition.boundary_conditions.vals[0].slope[1], 1);
}

TEST(EulerBeamHuangXMLTest, ValidatesSemanticsWithoutDtd)
{
  for (const auto& replacement : {
         std::pair<std::string, std::string>{"length=\"2\"", "length=\"-1\""},
         {"nodes=\"5\"", "nodes=\"1\""},
         {"Direction x=\"0\" y=\"2\"", "Direction x=\"0\" y=\"0\""},
         {"left=\"free\"", "left=\"unknown\""},
         {"anchorEnd=\"right\"", "anchorEnd=\"unknown\""}}) {
    try { parse(changed(replacement.first, replacement.second)); FAIL(); }
    catch (const XMLError& error) {
      EXPECT_EQ(error.status(), ELFF_XML_INVALID_VALUE);
      EXPECT_EQ(error.diagnostic().line, 1u);
    }
  }
}

TEST(EulerBeamHuangXMLTest, ReportsAttributeLocations)
{
  try { parse(changed("nodes=\"5\"", "nodes=\"bad\"")); FAIL(); }
  catch (const XMLError& error) {
    EXPECT_EQ(error.status(), ELFF_XML_INVALID_VALUE);
    EXPECT_EQ(error.diagnostic().line, 2u);
  }
}
#else
TEST(EulerBeamHuangXMLTest, ReportsUnavailable)
{
  // Disabled C stubs require a non-null handle but never dereference it.
  const char sentinel = 0;
  const auto node = reinterpret_cast<elff_xml_node_t>(&sentinel);
  try { parse_euler_beam_huang_xml(node); FAIL(); }
  catch (const XMLError& error) { EXPECT_EQ(error.status(), ELFF_XML_UNAVAILABLE); }
}
#endif
