#include "elff/c/models/ibm/IBCoupling.h"
#include <gtest/gtest.h>

#include "elff/c/models/beam/IBEulerBeamHuang.h"
#include "elff/c/models/ibm/IBForceCoupled.h"
#include "elff/c/models/ModelRuntime.h"
#include "elff/config/_config.hpp"

#include "elff/io/C/xmlCase.h"
#include "elff/models/rigidbody/IBRigidBodyFibre.hpp"
#include <cmath>
#include <vector>
#include <filesystem>
#include <fstream>
#include <string>

namespace {

class TemporaryCase
{
public:
  TemporaryCase(const char* name, const std::string& contents)
    : path_(name)
  {
    std::ofstream output(path_, std::ios::binary);
    output << contents;
  }

  ~TemporaryCase()
  {
    std::error_code error;
    std::filesystem::remove(path_, error);
  }

  const char* path() const { return path_.c_str(); }

private:
  std::filesystem::path path_;
};

std::string huang(const char* nodes,
                  const char* length = "1",
                  const char* direction_x = "1")
{
  return std::string(R"XML(
    <EulerBeamHuang>
      <Geometry length=")XML") + length + R"XML(" nodes=")XML" + nodes +
         R"XML("/>
      <Material bendingRigidity="0.001" linearDensity="1.5"/>
      <InitialLine anchorEnd="right">
        <Anchor x="0" y="0"/>
        <Direction x=")XML" + direction_x + R"XML(" y="0"/>
      </InitialLine>
      <BoundaryConditions left="free" right="simple"/>
    </EulerBeamHuang>
  )XML";
}

std::string case_xml(const std::string& models)
{
  return R"XML(<?xml version="1.0"?>
<Case formatVersion="1">
  <ELFF schemaVersion="1" dimensions="2" units="nondimensional">
    <Models>
)XML" + models + R"XML(
    </Models>
  </ELFF>
  <Basilisk><IgnoredByELFF/></Basilisk>
</Case>)XML";
}

#ifdef ELFF_USE_LIBXML2

TEST(ModelRuntimeXMLTest, BuildsRepeatedHuangModelsInDocumentOrder)
{
  TemporaryCase file("elff_runtime_ordered_models.xml",
                     case_xml(std::string("<!-- first -->") + huang("4") +
                              "\n<!-- second -->\n" + huang("7")));
  model_runtime_t runtime = nullptr;
  model_runtime_error error = {};
  ASSERT_EQ(model_runtime_create_from_case_file(file.path(), &runtime, &error), 0)
    << error.message;
  ASSERT_NE(runtime, nullptr);

  ib_model_view_array_t views = {};
  ASSERT_EQ(ib_coupling_get_model_views(ib_coupling_attach(runtime), &views), 0);
  ASSERT_EQ(views.count, 2U);
  EXPECT_EQ(views.data[0].kind, IB_MODEL_VIEW_FORCE_COUPLED);
  EXPECT_EQ(views.data[1].kind, IB_MODEL_VIEW_FORCE_COUPLED);
  EXPECT_EQ(ib_force_coupled_get_number_of_nodes(views.data[0].handle), 4);
  EXPECT_EQ(ib_force_coupled_get_number_of_nodes(views.data[1].handle), 7);

  model_runtime_delete(runtime);
}

TEST(ModelRuntimeXMLTest, UsesEmbeddedDtdAndRejectsLegacyModelNames)
{
  TemporaryCase file(
    "elff_runtime_legacy_model.xml",
    case_xml(std::string("<Model id=\"named\">") + huang("4") +
             "</Model>"));
  model_runtime_t runtime = reinterpret_cast<model_runtime_t>(1);
  model_runtime_error error = {};
  EXPECT_NE(model_runtime_create_from_case_file(file.path(), &runtime, &error), 0);
  EXPECT_EQ(runtime, nullptr);
  EXPECT_NE(error.message[0], '\0');
}

TEST(ModelRuntimeXMLTest, ReportsIndexedSemanticFailureTransactionally)
{
  TemporaryCase file("elff_runtime_invalid_second_model.xml",
                     case_xml(huang("4") + huang("7", "1", "0")));
  model_runtime_t runtime = reinterpret_cast<model_runtime_t>(1);
  model_runtime_error error = {};
  EXPECT_NE(model_runtime_create_from_case_file(file.path(), &runtime, &error), 0);
  EXPECT_EQ(runtime, nullptr);
  EXPECT_EQ(error.model_index, 1U);
  EXPECT_GT(error.line, 0);
  EXPECT_NE(std::string(error.message).find("direction"), std::string::npos);
}

TEST(ModelRuntimeXMLTest, RejectsNonpositiveGeometryWithModelLocation)
{
  TemporaryCase file("elff_runtime_invalid_geometry.xml",
                     case_xml(huang("4", "-1")));
  model_runtime_t runtime = nullptr;
  model_runtime_error error = {};
  EXPECT_NE(model_runtime_create_from_case_file(file.path(), &runtime, &error), 0);
  EXPECT_EQ(runtime, nullptr);
  EXPECT_EQ(error.model_index, 0U);
  EXPECT_GT(error.line, 0);
}

TEST(ModelRuntimeXMLTest, RejectsNonfiniteGeometryAndAnchorValues)
{
  TemporaryCase infinite_length("elff_runtime_infinite_length.xml",
                                case_xml(huang("4", "inf")));
  model_runtime_t runtime = nullptr;
  model_runtime_error error = {};
  EXPECT_NE(model_runtime_create_from_case_file(
              infinite_length.path(), &runtime, &error),
            0);
  EXPECT_EQ(runtime, nullptr);
  EXPECT_EQ(error.model_index, 0U);
  EXPECT_GT(error.line, 0);

  std::string model = huang("4");
  const std::string finite_anchor = "<Anchor x=\"0\" y=\"0\"/>";
  model.replace(model.find(finite_anchor),
                finite_anchor.size(),
                "<Anchor x=\"nan\" y=\"0\"/>");
  TemporaryCase nan_anchor("elff_runtime_nan_anchor.xml", case_xml(model));
  error = {};
  EXPECT_NE(model_runtime_create_from_case_file(
              nan_anchor.path(), &runtime, &error),
            0);
  EXPECT_EQ(runtime, nullptr);
  EXPECT_EQ(error.model_index, 0U);
  EXPECT_GT(error.line, 0);
}

TEST(ModelRuntimeXMLTest, ReturnedRuntimeIsAlreadySealed)
{
  TemporaryCase file("elff_runtime_sealed.xml", case_xml(huang("4")));
  model_runtime_t runtime = nullptr;
  model_runtime_error error = {};
  ASSERT_EQ(model_runtime_create_from_case_file(file.path(), &runtime, &error), 0)
    << error.message;

  ib_euler_beam_bcs_t bcs = {};
  bcs.end[0] = IB_EULER_BEAM_BC_LEFT;
  bcs.end[1] = IB_EULER_BEAM_BC_RIGHT;
  bcs.type[0] = IB_EULER_BEAM_BC_FREE;
  bcs.type[1] = IB_EULER_BEAM_BC_FREE;
  const vertex_t origin = {0, 0, 0};
  EXPECT_EQ(ib_euler_beam_huang_new(
              ib_coupling_attach(runtime), origin, bcs, 1.0, 1.0, 1.0, 4),
            nullptr);
  model_runtime_delete(runtime);
}


class RigidFibreXMLTest : public ::testing::Test {
protected:
  model_runtime_t runtime = nullptr;
  model_runtime_error error{};
  ~RigidFibreXMLTest() override { model_runtime_delete(runtime); }
  int load(const std::string& models, const std::string& domain = "", int dimension = 3) {
    const std::string xml = "<Case formatVersion='2' dimension='" + std::to_string(dimension) +
      "'><ELFF>" + domain + "<Models distribution='round-robin'>" + models + "</Models></ELFF></Case>";
    elff_xml_document_t document = nullptr;
    elff_xml_error diagnostic{};
    if (elff_xml_document_load_memory(xml.data(), xml.size(), &document, &diagnostic) != ELFF_XML_OK) return -1;
    const int result = model_runtime_create_from_case_document(document, &runtime, &error);
    elff_xml_document_destroy(document);
    return result;
  }
  std::string fibre() {
    return "<RigidBodyFibre><RigidFibreGeometry length='2' diameter='0.2' nodes='5'/>"
      "<RigidFibreMaterial><LinearDensity>3</LinearDensity></RigidFibreMaterial>"
      "<RigidBodyInitialCondition><Center x='0' y='0' z='0'/>"
      "<Orientation w='2' x='0' y='0' z='2'/><Velocity x='1' y='2' z='3'/>"
      "<AngularVelocity x='0' y='0' z='1'/></RigidBodyInitialCondition></RigidBodyFibre>";
  }
};

TEST_F(RigidFibreXMLTest, ConstructsFibreWithOrientationVelocityAndMeasures) {
  ASSERT_EQ(load(fibre()), 0) << error.message;
  ib_model_view_array_t views{};
  ASSERT_EQ(ib_coupling_get_model_views(ib_coupling_attach(runtime), &views), 0);
  ASSERT_EQ(views.count, 1u);
  auto* fibre = dynamic_cast<ELFF::Models::IBRigidBodyFibre*>(
      static_cast<ELFF::Models::IBForceCoupled*>(views.data[0].handle));
  ASSERT_NE(fibre, nullptr);
  EXPECT_NEAR(fibre->total_mass(), 6, 1e-12);
  auto mesh = ib_force_coupled_get_current(views.data[0].handle);
  ASSERT_EQ(mesh.n, 5);
  ELFF::real_t measure = 0;
  for (int i = 0; i < mesh.n; ++i) measure += mesh.measure[i];
  EXPECT_NEAR(measure, 2, 1e-12);
  EXPECT_NEAR(mesh.position[0].x, 0, 1e-12);
  EXPECT_NEAR(mesh.position[0].y, -1, 1e-12);
  EXPECT_NEAR(mesh.velocity[0].x, 2, 1e-12);
  EXPECT_NEAR(mesh.velocity[0].y, 2, 1e-12);
  EXPECT_NEAR(mesh.velocity[0].z, 3, 1e-12);
  ib_mesh_free(&mesh);
}

TEST_F(RigidFibreXMLTest, RejectsInvalidInputsBeforeConstruction) {
  for (const auto& change : std::vector<std::pair<std::string, std::string>>{
      {"length='2'", "length='0'"}, {"length='2'", "length='1e200'"},
      {"diameter='0.2'", "diameter='0'"}, {"diameter='0.2'", "diameter='1e-10'"},
      {"nodes='5'", "nodes='1'"}, {"<LinearDensity>3", "<LinearDensity>-1"},
      {"x='0'", "x='nan'"}}) {
    auto xml = fibre();
    xml.replace(xml.find(change.first), change.first.size(), change.second);
    EXPECT_NE(load(xml), 0) << change.second;
    EXPECT_EQ(runtime, nullptr);
  }
  auto xml = fibre();
  xml.replace(xml.find("w='2'"), 5, "w='0'");
  xml.replace(xml.find("z='2'"), 5, "z='0'");
  EXPECT_NE(load(xml), 0);
  EXPECT_EQ(runtime, nullptr);
  EXPECT_NE(load(fibre(), "", 2), 0);
  EXPECT_EQ(runtime, nullptr);
}

TEST_F(RigidFibreXMLTest, InsertsRepeatedRigidFibres) {
  auto xml = fibre();
  const std::string center = "<Center x='0' y='0' z='0'/>";
  xml.replace(xml.find(center), center.size(),
      "<RandomCenter><Min x='3' y='-3' z='-3'/><Max x='4' y='3' z='3'/></RandomCenter>");
  EXPECT_NE(load(xml), 0); // Random placement requires an insertion group.
  xml.insert(xml.find("</RigidBodyFibre>"), "<Insertion count='3' clearance='0.1' maxAttempts='100'/>");
  ASSERT_EQ(load(fibre() + xml,
      "<InsertionDomain><Min x='-5' y='-5' z='-5'/><Max x='5' y='5' z='5'/></InsertionDomain>"), 0)
      << error.message;
  ib_model_view_array_t views{};
  ASSERT_EQ(ib_coupling_get_model_views(ib_coupling_attach(runtime), &views), 0);
  ASSERT_EQ(views.count, 4u);
  EXPECT_EQ(ib_force_coupled_get_number_of_nodes(views.data[0].handle), 5);
  for (size_t i = 1; i < views.count; ++i)
    EXPECT_EQ(ib_force_coupled_get_number_of_nodes(views.data[i].handle), 5);
}

#else

TEST(ModelRuntimeXMLTest, ReportsUnavailableWithoutLibxml2)
{
  model_runtime_t runtime = nullptr;
  model_runtime_error error = {};
  EXPECT_NE(model_runtime_create_from_case_file("case.xml", &runtime, &error), 0);
  EXPECT_EQ(runtime, nullptr);
}

#endif

} // namespace
