#include <gtest/gtest.h>
#include "elff/config/_config.hpp"
#include "elff/io/case/models/EulerBeamHuangXML.hpp"
#include "elff/models/beam/EulerBeamHuangBuilder.hpp"
#include "elff/io/CXX/xmlCase.hpp"
#include <memory>
#include <string>

using namespace ELFF::Models;
using ELFF::IO::CXX::XMLError;
using ELFF::IO::CXX::XMLNode;
using ELFF::IO::Case::EulerBeamHuangXML;

#ifdef ELFF_USE_LIBXML2
namespace {
const std::string valid = R"(<EulerBeamHuang>
<Geometry length="2" nodes="5"/>
<Material bendingRigidity="0.01" linearDensity="1.5"/>
<InitialLine anchorEnd="right"><Anchor x="1" y="2"/><Direction x="0" y="2"/></InitialLine>
<BoundaryConditions left="free" right="clamped"/>
</EulerBeamHuang>)";

struct Parsed { EulerBeamHuangDefinition definition; ELFF::IO::Case::InitialLineTemplate placement; };
Parsed parse(const std::string& xml)
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
  const XMLNode node(root);
  auto definition = EulerBeamHuangXML::parse(node);
  return {definition, EulerBeamHuangXML::parse_placement(node, definition)};
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
  const auto parsed = parse(valid);
  const auto& definition = parsed.definition;
  auto beam = EulerBeamHuangBuilder::build(definition, parsed.placement.fixed);
  EXPECT_EQ(definition.nodes, 5u);
  EXPECT_EQ(definition.source_line, 1);
  EXPECT_DOUBLE_EQ(definition.length, 2);
  EXPECT_TRUE(parsed.placement.fixed.anchor_at_right);
  EXPECT_DOUBLE_EQ(beam->get_mesh().get_centerline(0)[1], 4);
  EXPECT_DOUBLE_EQ(beam->get_mesh().get_centerline(4)[1], 2);
  EXPECT_DOUBLE_EQ(beam->get_mesh().get_slope(0)[1], -1);
  EXPECT_EQ(definition.boundary_conditions.type[1], EulerBeam::clamped_bc);
}

TEST(EulerBeamHuangXMLTest, SupportsLeftAnchoring)
{
  const auto parsed = parse(changed("anchorEnd=\"right\"", "anchorEnd=\"left\""));
  auto beam = EulerBeamHuangBuilder::build(parsed.definition, parsed.placement.fixed);
  EXPECT_FALSE(parsed.placement.fixed.anchor_at_right);
  EXPECT_DOUBLE_EQ(beam->get_mesh().get_centerline(0)[1], 2);
  EXPECT_DOUBLE_EQ(beam->get_mesh().get_centerline(4)[1], 4);
  EXPECT_DOUBLE_EQ(beam->get_mesh().get_slope(0)[1], 1);
}

TEST(EulerBeamHuangXMLTest, ValidatesSemanticsWithoutDtd)
{
  const unsigned lines[] = {1, 2, 4, 5, 4};
  size_t index = 0;
  for (const auto& replacement : {
         std::pair<std::string, std::string>{"length=\"2\"", "length=\"-1\""},
         {"nodes=\"5\"", "nodes=\"1\""},
         {"Direction x=\"0\" y=\"2\"", "Direction x=\"0\" y=\"0\""},
         {"left=\"free\"", "left=\"unknown\""},
         {"anchorEnd=\"right\"", "anchorEnd=\"unknown\""}}) {
    try { parse(changed(replacement.first, replacement.second)); FAIL(); }
    catch (const XMLError& error) {
      EXPECT_EQ(error.status(), ELFF_XML_INVALID_VALUE);
      EXPECT_EQ(error.diagnostic().line, lines[index]);
    }
    ++index;
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
  try { EulerBeamHuangXML::parse(XMLNode(node)); FAIL(); }
  catch (const XMLError& error) { EXPECT_EQ(error.status(), ELFF_XML_UNAVAILABLE); }
}
#endif

TEST(EulerBeamHuangBuilderTest, BuildsWithoutXml)
{
  EulerBeamHuangDefinition definition;
  definition.length = 2;
  definition.nodes = 5;
  definition.bending_rigidity = 0.01;
  definition.linear_density = 1.5;
  definition.boundary_conditions.type[0] = EulerBeam::free_bc;
  definition.boundary_conditions.type[1] = EulerBeam::clamped_bc;
  ELFF::IO::Case::InitialLineDefinition placement;
  placement.anchor = {1, 2, 3};
  placement.direction = {0, 2, 0};
  ASSERT_EQ(EulerBeamHuangBuilder::validation_error(definition, placement), nullptr);
  auto beam = EulerBeamHuangBuilder::build(definition, placement);
  EXPECT_DOUBLE_EQ(beam->get_mesh().get_centerline(4)[1], 4);
  EXPECT_DOUBLE_EQ(beam->get_mesh().get_slope(0)[1], 1);
  placement.direction = {0, 0, 0};
  EXPECT_NE(EulerBeamHuangBuilder::validation_error(definition, placement), nullptr);
}
