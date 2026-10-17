#include <gtest/gtest.h>
#include <elff/config/_config.hpp>
#include "basilisk/library/io/input-case.h"
#include <filesystem>
#include <fstream>
#include <string>

#ifdef ELFF_USE_LIBXML2
namespace {
const std::string settings = R"xml(
<Basilisk caseType="unrelated_solver">
  <Domain size="8" originX="-4" originY="-4"/>
  <Time step="0.001" end="1"/>
  <Fluid density="1" dynamicViscosity="0.01"><InletVelocity x="1" y="0"/></Fluid>
  <Grid minLevel="3" maxLevel="6"/>
  <Output directory="output" interval="0.1"/>
  <UnrelatedExtension/>
</Basilisk>)xml";

class BasiliskCaseReaderTest : public ::testing::Test {
protected:
  elff_xml_document_t document = nullptr;
  elff_xml_error error{};
  BasiliskCaseConfiguration config{};
  ~BasiliskCaseReaderTest() override { elff_xml_document_destroy(document); }
  elff_xml_status read(const std::string& xml) {
    elff_xml_document_destroy(document);
    document = nullptr;
    auto status = elff_xml_document_load_memory(xml.data(), xml.size(), &document, &error);
    return status == ELFF_XML_OK ? basilisk_case_read(document, &config, &error) : status;
  }
};

TEST_F(BasiliskCaseReaderTest, ReadsCommonSettingsWithoutAnyFibreOrElffModel) {
  ASSERT_EQ(read("<Case>" + settings + "</Case>"), ELFF_XML_OK) << error.message;
  EXPECT_DOUBLE_EQ(config.domain_size, 8);
  EXPECT_DOUBLE_EQ(config.origin[2], 0);
  EXPECT_DOUBLE_EQ(config.dynamic_viscosity, .01);
  EXPECT_DOUBLE_EQ(config.inlet_velocity[0], 1);
  EXPECT_EQ(config.immersed_boundary_level, 6);
  EXPECT_EQ(config.has_uniform_force_density, 0);
  EXPECT_DOUBLE_EQ(config.uniform_force_density[0], 0);
  EXPECT_STREQ(config.output_directory, "output");
}

TEST_F(BasiliskCaseReaderTest, ReadsOptionalImmersedBoundarySettings) {
  std::string xml = settings;
  xml.insert(xml.find("</Basilisk>"),
    "<ImmersedBoundary level='7'><UniformForceDensity x='0' y='-2'/></ImmersedBoundary>");
  ASSERT_EQ(read("<Case>" + xml + "</Case>"), ELFF_XML_OK) << error.message;
  EXPECT_EQ(config.immersed_boundary_level, 7);
  EXPECT_EQ(config.has_uniform_force_density, 1);
  EXPECT_DOUBLE_EQ(config.uniform_force_density[1], -2);
  EXPECT_DOUBLE_EQ(config.uniform_force_density[2], 0);
}

TEST_F(BasiliskCaseReaderTest, RejectsInvalidCommonSettings) {
  for (const auto& replacement : {
      std::pair<std::string, std::string>{"size=\"8\"", "size=\"-1\""},
      {"step=\"0.001\"", "step=\"0\""},
      {"density=\"1\"", "density=\"nan\""},
      {"maxLevel=\"6\"", "maxLevel=\"2\""}}) {
    auto xml = settings;
    xml.replace(xml.find(replacement.first), replacement.first.size(), replacement.second);
    EXPECT_NE(read("<Case>" + xml + "</Case>"), ELFF_XML_OK) << replacement.second;
  }
}

TEST_F(BasiliskCaseReaderTest, LoadReturnsRenderedDocumentForAdditionalReaders) {
  const std::filesystem::path path = "basilisk_case_reader_load.xml";
  struct Cleanup {
    std::filesystem::path path;
    ~Cleanup() { std::filesystem::remove(path); }
  } cleanup{path};
  auto xml = settings;
  xml.replace(xml.find("size=\"8\""), 8, "size=\"${L}\"");
  {
    std::ofstream file(path);
    file << "<Case><Parameters><Parameter name='L' value='8'/></Parameters>"
         << xml << "<Other length='${L/2}'/></Case>";
  }
  ASSERT_EQ(basilisk_case_load(path.c_str(), 0, &document, &config, &error), ELFF_XML_OK)
    << error.message;
  EXPECT_DOUBLE_EQ(config.domain_size, 8);
  elff_xml_node_t root = nullptr, other = nullptr;
  ASSERT_EQ(elff_xml_document_root(document, "Case", &root, &error), ELFF_XML_OK);
  ASSERT_EQ(elff_xml_node_unique_child(root, "Other", &other, &error), ELFF_XML_OK);
  double length = 0;
  ASSERT_EQ(elff_xml_node_attribute_double(other, "length", &length, &error), ELFF_XML_OK);
  EXPECT_DOUBLE_EQ(length, 4);
  elff_xml_document_destroy(document);
  document = nullptr;
  { std::ofstream file(path); file << "<Case><Basilisk/></Case>"; }
  EXPECT_NE(basilisk_case_load(path.c_str(), 0, &document, &config, &error), ELFF_XML_OK);
  EXPECT_EQ(document, nullptr);
}
} // namespace
#endif
