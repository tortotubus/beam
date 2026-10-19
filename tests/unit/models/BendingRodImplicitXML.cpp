#include <gtest/gtest.h>
#include "elff/models/articulated/BendingRodImplicitBuilder.hpp"
#include "elff/models/ibm/IBRuntime.hpp"
#include "elff/c/models/ibm/IBForceCoupled.h"
#include "elff/io/CXX/xmlCase.hpp"
#include "embeddedELFFDtd.hpp"
#include <filesystem>
#include <fstream>
#include <memory>
#include <string>

using namespace ELFF::Models;
using ELFF::IO::CXX::XMLError;
#ifdef ELFF_USE_LIBXML2
namespace {
const std::string valid = R"(<BendingRodSegmentsImplicit>
<Geometry length="2" nodes="5"/>
<Material bendingRigidity="0.01" linearDensity="1.5"/>
<RodInitialLine><Anchor x="1" y="2"/><Direction x="0" y="2"/></RodInitialLine>
<RodBoundary type="pinned"/>
</BendingRodSegmentsImplicit>)";
using Document = std::unique_ptr<elff_xml_document_handle, decltype(&elff_xml_document_destroy)>;
Document document(const std::string &xml) {
  elff_xml_document_t raw = nullptr;
  elff_xml_error error{};
  if (elff_xml_document_load_memory(xml.data(), xml.size(), &raw, &error) != ELFF_XML_OK)
    throw std::runtime_error(error.message);
  return Document(raw, elff_xml_document_destroy);
}
elff_xml_node_t root(elff_xml_document_t doc, const char *name) {
  elff_xml_node_t node = nullptr;
  elff_xml_error error{};
  if (elff_xml_document_root(doc, name, &node, &error) != ELFF_XML_OK)
    throw std::runtime_error(error.message);
  return node;
}
BendingRodImplicitDefinition parse(const std::string &xml) {
  auto doc = document(xml);
  return parse_bending_rod_implicit_xml(root(doc.get(), "BendingRodSegmentsImplicit"));
}
std::string changed(const std::string &from, const std::string &to) {
  auto xml = valid;
  xml.replace(xml.find(from), from.size(), to);
  return xml;
}
std::string section(const std::string &models) {
  return "<ELFF schemaVersion='1' dimensions='3' units='SI'><Models>" + models + "</Models></ELFF>";
}
std::string case_xml(const std::string &models) {
  return "<Case formatVersion='1'>" + section(models) + "</Case>";
}
using Runtime = std::unique_ptr<void, decltype(&ib_runtime_delete)>;
}

TEST(BendingRodImplicitXMLTest, DefaultsAndGeometrySurviveDocumentDestruction)
{
  for (const auto &boundary : {"pinned", "free"}) {
    const auto definition = parse(changed("type=\"pinned\"", std::string("type=\"") + boundary + "\""));
    EXPECT_EQ(definition.config.links, 4);
    EXPECT_DOUBLE_EQ(definition.config.total_mass, 3);
    EXPECT_DOUBLE_EQ(definition.config.spectral_radius, 0.8);
    EXPECT_DOUBLE_EQ(definition.config.nonlinear_tolerance, 1e-9);
    EXPECT_EQ(definition.config.max_iterations, 20);
    EXPECT_DOUBLE_EQ(definition.config.gravity.y(), -9.81);
    auto rod = build_bending_rod_implicit(definition);
    const auto &mesh = rod->GetCurrent();
    EXPECT_DOUBLE_EQ(mesh.GetPoints().front().x, 1);
    EXPECT_DOUBLE_EQ(mesh.GetPoints().front().y, 2);
    EXPECT_DOUBLE_EQ(mesh.GetPoints().back().y, 4);
    EXPECT_DOUBLE_EQ(mesh.GetMeasures().front(), 0.25);
    EXPECT_DOUBLE_EQ(mesh.GetMeasures()[1], 0.5);
    EXPECT_DOUBLE_EQ(mesh.GetVelocity().back().y, 0);
  }
}

TEST(BendingRodImplicitXMLTest, ControlsAndSemanticFailuresWithoutDtd)
{
  auto xml = changed("</BendingRodSegmentsImplicit>",
    "<Gravity x='1' y='2' z='3'/><ImplicitSolver spectralRadius='0.5' nonlinearTolerance='1e-12' maxIterations='30'/></BendingRodSegmentsImplicit>");
  auto definition = parse(xml);
  EXPECT_DOUBLE_EQ(definition.config.gravity.z(), 3);
  EXPECT_DOUBLE_EQ(definition.config.spectral_radius, 0.5);
  EXPECT_DOUBLE_EQ(definition.config.nonlinear_tolerance, 1e-12);
  EXPECT_EQ(definition.config.max_iterations, 30);
  EXPECT_EQ(parse(changed("</BendingRodSegmentsImplicit>", "<ImplicitSolver/></BendingRodSegmentsImplicit>")).config.max_iterations, 20);
  for (const auto &change : {
      std::pair<std::string, std::string>{"nodes=\"5\"", "nodes=\"1\""},
      {"nodes=\"5\"", "nodes=\"2147483647\""},
      {"length=\"2\"", "length=\"-1\""},
      {"length=\"2\"", "length=\"1e300\""},
      {"linearDensity=\"1.5\"", "linearDensity=\"0\""},
      {"bendingRigidity=\"0.01\"", "bendingRigidity=\"-1\""},
      {"Direction x=\"0\" y=\"2\"", "Direction x=\"0\" y=\"0\""},
      {"type=\"pinned\"", "type=\"clamped\""},
      {"</BendingRodSegmentsImplicit>", "<ImplicitSolver spectralRadius='2'/></BendingRodSegmentsImplicit>"},
      {"</BendingRodSegmentsImplicit>", "<ImplicitSolver nonlinearTolerance='0'/></BendingRodSegmentsImplicit>"},
      {"</BendingRodSegmentsImplicit>", "<ImplicitSolver maxIterations='0'/></BendingRodSegmentsImplicit>"},
      {"</BendingRodSegmentsImplicit>", "<Gravity x='NaN' y='0'/></BendingRodSegmentsImplicit>"},
      {"</BendingRodSegmentsImplicit>", "<ImplicitSolver/><ImplicitSolver/></BendingRodSegmentsImplicit>"}}) {
    SCOPED_TRACE(change.second);
    EXPECT_THROW(parse(changed(change.first, change.second)), XMLError);
  }
  EXPECT_NO_THROW(parse(changed("bendingRigidity=\"0.01\"", "bendingRigidity=\"0\"")));
  try { parse(changed("nodes=\"5\"", "nodes=\"bad\"")); FAIL(); }
  catch (const XMLError &error) { EXPECT_EQ(error.diagnostic().line, 2u); }
}

TEST(BendingRodImplicitXMLTest, EmbeddedDtdAcceptsMixedModelsAndRejectsBadStructure)
{
  const std::string huang = R"(<EulerBeamHuang><Geometry length='1' nodes='3'/>
<Material bendingRigidity='1' linearDensity='1'/><InitialLine anchorEnd='left'>
<Anchor x='0' y='0'/><Direction x='1' y='0'/></InitialLine>
<BoundaryConditions left='simple' right='free'/></EulerBeamHuang>)";
  for (bool invalid : {false, true}) {
    auto doc = document(section(huang + (invalid ? changed("RodBoundary", "BadBoundary") : valid)));
    char error[512]{};
    const auto rc = elff_xml_dtd_validate_node(root(doc.get(), "ELFF"),
      reinterpret_cast<const char*>(elff_embedded_elff_v1_dtd), elff_embedded_elff_v1_dtd_size,
      error, sizeof(error));
    EXPECT_EQ(rc, invalid ? ELFF_XML_DTD_INVALID_DOCUMENT : ELFF_XML_DTD_VALID) << error;
  }
}

TEST(BendingRodImplicitXMLTest, RuntimeFileDocumentSteppingAndCheckpoint)
{
  const auto xml = changed("</BendingRodSegmentsImplicit>", "<Gravity x='0' y='0'/></BendingRodSegmentsImplicit>");
  auto case_text = case_xml(xml);
  case_text.insert(case_text.find('>') + 1, "<Parameters><Parameter name='nodes' value='5'/></Parameters>");
  case_text.replace(case_text.find("nodes=\"5\""), 9, "nodes=\"${parameters.nodes}\"");
  const auto file = std::filesystem::temp_directory_path() / "elff_implicit_rod_case.xml";
  { std::ofstream output(file); output << case_text; }
  auto source = document(case_text);
  elff_xml_document_t rendered = nullptr;
  elff_xml_error xml_error{};
  ASSERT_EQ(elff_xml_case_render(source.get(), &rendered, &xml_error), ELFF_XML_OK) << xml_error.message;
  Document doc(rendered, elff_xml_document_destroy);
  source.reset();
  ib_runtime_t raw = nullptr, raw_restored = nullptr;
  ib_runtime_error error{};
  ASSERT_EQ(ib_runtime_create_from_case_file(file.c_str(), &raw, &error), 0) << error.message;
  Runtime runtime(raw, ib_runtime_delete);
  std::filesystem::remove(file);
  ASSERT_EQ(ib_runtime_create_from_case_document(doc.get(), &raw_restored, &error), 0) << error.message;
  Runtime restored(raw_restored, ib_runtime_delete);
  doc.reset();
  ib_model_view_array_t views{}, restored_views{};
  ASSERT_EQ(ib_runtime_get_model_views(raw, &views), 0);
  ASSERT_EQ(ib_runtime_get_model_views(raw_restored, &restored_views), 0);
  ASSERT_EQ(views.count, 1u);
  EXPECT_EQ(views.data[0].kind, IB_MODEL_VIEW_FORCE_COUPLED);
  auto direct = build_bending_rod_implicit(parse(xml));
  std::vector<vertex_t> force(5, {0.2, -0.1, 0.3});
  std::vector<IBMesh::IBVertex> cpp_force(5, {0.2, -0.1, 0.3});
  char diagnostic[512]{};
  for (int step = 0; step < 3; ++step) {
    ib_mesh_t next{};
    ASSERT_EQ(ib_force_coupled_get_next_checked(views.data[0].handle, force.data(), 5, 0.01,
      &next, diagnostic, sizeof(diagnostic)), 0) << diagnostic;
    auto &expected = direct->GetNext(cpp_force, 0.01);
    EXPECT_NEAR(next.position[4].x, expected.GetPoints()[4].x, 1e-13);
    EXPECT_NEAR(next.velocity[4].z, expected.GetVelocity()[4].z, 1e-13);
    ib_mesh_free(&next);
  }
  const auto checkpoint = std::filesystem::temp_directory_path() / "elff_implicit_rod_checkpoint.bin";
  ASSERT_EQ(ib_runtime_checkpoint(raw, checkpoint.c_str()), 0);
  ASSERT_EQ(ib_runtime_restore(raw_restored, checkpoint.c_str()), 0);
  std::filesystem::remove(checkpoint);
  for (double dt : {0.01, 0.005}) {
    auto a = ib_force_coupled_get_next(views.data[0].handle, force.data(), 5, dt);
    auto b = ib_force_coupled_get_next(restored_views.data[0].handle, force.data(), 5, dt);
    ASSERT_EQ(a.n, 5); ASSERT_EQ(b.n, 5);
    EXPECT_DOUBLE_EQ(a.position[4].x, b.position[4].x);
    EXPECT_DOUBLE_EQ(a.velocity[4].z, b.velocity[4].z);
    ib_mesh_free(&a); ib_mesh_free(&b);
  }
}
TEST(BendingRodImplicitXMLTest, RuntimeReportsIndexedRodSemanticError)
{
  auto doc = document(case_xml(valid + changed("nodes=\"5\"", "nodes=\"1\"")));
  ib_runtime_t runtime = nullptr;
  ib_runtime_error error{};
  EXPECT_NE(ib_runtime_create_from_case_document(doc.get(), &runtime, &error), 0);
  EXPECT_EQ(runtime, nullptr);
  EXPECT_EQ(error.model_index, 1u);
  EXPECT_GT(error.line, 0);
  EXPECT_NE(error.message[0], '\0');
}
#else
TEST(BendingRodImplicitXMLTest, ReportsUnavailable)
{
  const char sentinel = 0;
  try { parse_bending_rod_implicit_xml(reinterpret_cast<elff_xml_node_t>(&sentinel)); FAIL(); }
  catch (const XMLError &error) { EXPECT_EQ(error.status(), ELFF_XML_UNAVAILABLE); }
}
#endif
