#include <gtest/gtest.h>

#include "elff/c/models/beam/IBEulerBeamHuang.h"
#include "elff/c/models/ibm/IBForceCoupled.h"
#include "elff/c/models/ibm/IBRuntime.h"
#include "elff/config/_config.hpp"

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

TEST(IBRuntimeXMLTest, BuildsRepeatedHuangModelsInDocumentOrder)
{
  TemporaryCase file("elff_runtime_ordered_models.xml",
                     case_xml(std::string("<!-- first -->") + huang("4") +
                              "\n<!-- second -->\n" + huang("7")));
  ib_runtime_t runtime = nullptr;
  ib_runtime_error error = {};
  ASSERT_EQ(ib_runtime_create_from_case_file(file.path(), &runtime, &error), 0)
    << error.message;
  ASSERT_NE(runtime, nullptr);

  ib_model_view_array_t views = {};
  ASSERT_EQ(ib_runtime_get_model_views(runtime, &views), 0);
  ASSERT_EQ(views.count, 2U);
  EXPECT_EQ(views.data[0].kind, IB_MODEL_VIEW_FORCE_COUPLED);
  EXPECT_EQ(views.data[1].kind, IB_MODEL_VIEW_FORCE_COUPLED);
  EXPECT_EQ(ib_force_coupled_get_number_of_nodes(views.data[0].handle), 4);
  EXPECT_EQ(ib_force_coupled_get_number_of_nodes(views.data[1].handle), 7);

  ib_runtime_delete(runtime);
}

TEST(IBRuntimeXMLTest, UsesEmbeddedDtdAndRejectsLegacyModelNames)
{
  TemporaryCase file(
    "elff_runtime_legacy_model.xml",
    case_xml(std::string("<Model id=\"named\">") + huang("4") +
             "</Model>"));
  ib_runtime_t runtime = reinterpret_cast<ib_runtime_t>(1);
  ib_runtime_error error = {};
  EXPECT_NE(ib_runtime_create_from_case_file(file.path(), &runtime, &error), 0);
  EXPECT_EQ(runtime, nullptr);
  EXPECT_NE(error.message[0], '\0');
}

TEST(IBRuntimeXMLTest, ReportsIndexedSemanticFailureTransactionally)
{
  TemporaryCase file("elff_runtime_invalid_second_model.xml",
                     case_xml(huang("4") + huang("7", "1", "0")));
  ib_runtime_t runtime = reinterpret_cast<ib_runtime_t>(1);
  ib_runtime_error error = {};
  EXPECT_NE(ib_runtime_create_from_case_file(file.path(), &runtime, &error), 0);
  EXPECT_EQ(runtime, nullptr);
  EXPECT_EQ(error.model_index, 1U);
  EXPECT_GT(error.line, 0);
  EXPECT_NE(std::string(error.message).find("direction"), std::string::npos);
}

TEST(IBRuntimeXMLTest, RejectsNonpositiveGeometryWithModelLocation)
{
  TemporaryCase file("elff_runtime_invalid_geometry.xml",
                     case_xml(huang("4", "-1")));
  ib_runtime_t runtime = nullptr;
  ib_runtime_error error = {};
  EXPECT_NE(ib_runtime_create_from_case_file(file.path(), &runtime, &error), 0);
  EXPECT_EQ(runtime, nullptr);
  EXPECT_EQ(error.model_index, 0U);
  EXPECT_GT(error.line, 0);
}

TEST(IBRuntimeXMLTest, RejectsNonfiniteGeometryAndAnchorValues)
{
  TemporaryCase infinite_length("elff_runtime_infinite_length.xml",
                                case_xml(huang("4", "inf")));
  ib_runtime_t runtime = nullptr;
  ib_runtime_error error = {};
  EXPECT_NE(ib_runtime_create_from_case_file(
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
  EXPECT_NE(ib_runtime_create_from_case_file(
              nan_anchor.path(), &runtime, &error),
            0);
  EXPECT_EQ(runtime, nullptr);
  EXPECT_EQ(error.model_index, 0U);
  EXPECT_GT(error.line, 0);
}

TEST(IBRuntimeXMLTest, ReturnedRuntimeIsAlreadySealed)
{
  TemporaryCase file("elff_runtime_sealed.xml", case_xml(huang("4")));
  ib_runtime_t runtime = nullptr;
  ib_runtime_error error = {};
  ASSERT_EQ(ib_runtime_create_from_case_file(file.path(), &runtime, &error), 0)
    << error.message;

  ib_euler_beam_bcs_t bcs = {};
  bcs.end[0] = IB_EULER_BEAM_BC_LEFT;
  bcs.end[1] = IB_EULER_BEAM_BC_RIGHT;
  bcs.type[0] = IB_EULER_BEAM_BC_FREE;
  bcs.type[1] = IB_EULER_BEAM_BC_FREE;
  const vertex_t origin = {0, 0, 0};
  EXPECT_EQ(ib_euler_beam_huang_new(
              runtime, origin, bcs, 1.0, 1.0, 1.0, 4),
            nullptr);
  ib_runtime_delete(runtime);
}

#else

TEST(IBRuntimeXMLTest, ReportsUnavailableWithoutLibxml2)
{
  ib_runtime_t runtime = nullptr;
  ib_runtime_error error = {};
  EXPECT_NE(ib_runtime_create_from_case_file("case.xml", &runtime, &error), 0);
  EXPECT_EQ(runtime, nullptr);
}

#endif

} // namespace
