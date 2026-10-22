#include <gtest/gtest.h>

#include <elff/config/_config.hpp>
#include <elff/io/C/xmlCase.h>
#include "embeddedELFFDtd.hpp"

#include <fstream>
#include <iterator>
#include <string>
#include <vector>

namespace {

#ifdef ELFF_USE_LIBXML2

TEST(XMLNumericTextTest, ReadsLeafNumbersAndPreservesOutputOnFailure)
{
  for (const std::string content : {" 1.25 ", "", "NaN", "1 trailing", "<Nested>1</Nested>", "${x}"}) {
    const std::string xml = "<Value>" + content + "</Value>";
    elff_xml_document_t doc = nullptr;
    elff_xml_node_t node = nullptr;
    elff_xml_error error{};
    ASSERT_EQ(elff_xml_document_load_memory(xml.data(), xml.size(), &doc, &error), ELFF_XML_OK);
    ASSERT_EQ(elff_xml_document_root(doc, "Value", &node, &error), ELFF_XML_OK);
    double value = 42;
    const auto status = elff_xml_node_text_double(node, &value, &error);
    if (content == " 1.25 ") {
      EXPECT_EQ(status, ELFF_XML_OK);
      EXPECT_DOUBLE_EQ(value, 1.25);
    } else {
      EXPECT_EQ(status, ELFF_XML_INVALID_VALUE);
      EXPECT_DOUBLE_EQ(value, 42);
      EXPECT_GT(error.line, 0u);
    }
    elff_xml_document_destroy(doc);
  }
}

class XMLCaseDocumentTest : public ::testing::Test
{
protected:
  void SetUp() override
  {
    elff_xml_error error = {};
    ASSERT_EQ(elff_xml_case_load_rendered_file(
                ELFF_FLAPPING_FIBRE_CASE_PATH, 0, &document_, &error),
              ELFF_XML_OK)
      << error.message;
  }

  void TearDown() override { elff_xml_document_destroy(document_); }

  elff_xml_document_t document_ = nullptr;
};

TEST_F(XMLCaseDocumentTest, ReadsTypedValuesFromCanonicalCase)
{
  elff_xml_error error = {};
  elff_xml_node_t root = nullptr;
  elff_xml_node_t basilisk = nullptr;
  elff_xml_node_t domain = nullptr;
  elff_xml_node_t grid = nullptr;
  int format_version = 0;
  int min_level = 0;
  double domain_size = 0.0;

  ASSERT_EQ(elff_xml_document_root(document_, "Case", &root, &error),
            ELFF_XML_OK)
    << error.message;
  EXPECT_EQ(elff_xml_node_attribute_int(
              root, "formatVersion", &format_version, &error),
            ELFF_XML_OK)
    << error.message;
  EXPECT_EQ(format_version, 2);

  ASSERT_EQ(elff_xml_node_unique_child(
              root, "Basilisk", &basilisk, &error),
            ELFF_XML_OK)
    << error.message;
  ASSERT_EQ(elff_xml_node_unique_child(basilisk, "Grid", &grid, &error),
            ELFF_XML_OK)
    << error.message;

  ASSERT_EQ(elff_xml_node_unique_child(grid, "Length", &domain, &error), ELFF_XML_OK);
  EXPECT_EQ(elff_xml_node_attribute_double(
              domain, "x", &domain_size, &error),
            ELFF_XML_OK)
    << error.message;
  EXPECT_DOUBLE_EQ(domain_size, 16.0);
  elff_xml_node_t adaptation = nullptr;
  ASSERT_EQ(elff_xml_node_unique_child(grid, "AdaptWavelet", &adaptation, &error), ELFF_XML_OK);
  EXPECT_EQ(elff_xml_node_attribute_int(
              adaptation, "minLevel", &min_level, &error),
            ELFF_XML_OK)
    << error.message;
  EXPECT_EQ(min_level, 6);
}

TEST_F(XMLCaseDocumentTest, ValidatesV2CaseWithEmbeddedSchema)
{
  elff_xml_error error{};
  EXPECT_EQ(elff_xml_case_validate(document_, &error), ELFF_XML_OK) << error.message;
}

TEST_F(XMLCaseDocumentTest, ReportsRequiredStringCapacity)
{
  elff_xml_error error = {};
  elff_xml_node_t root = nullptr;
  elff_xml_node_t basilisk = nullptr;
  elff_xml_node_t output = nullptr;
  size_t required = 0;

  ASSERT_EQ(elff_xml_document_root(document_, "Case", &root, &error),
            ELFF_XML_OK);
  ASSERT_EQ(elff_xml_node_unique_child(
              root, "Basilisk", &basilisk, &error),
            ELFF_XML_OK);
  ASSERT_EQ(elff_xml_node_unique_child(
              basilisk, "Output", &output, &error),
            ELFF_XML_OK);

  EXPECT_EQ(elff_xml_node_attribute_string(
              output, "directory", nullptr, 0, &required, &error),
            ELFF_XML_BUFFER_TOO_SMALL);
  ASSERT_GT(required, 1U);

  std::vector<char> value(required);
  EXPECT_EQ(elff_xml_node_attribute_string(output,
                                            "directory",
                                            value.data(),
                                            value.size(),
                                            &required,
                                            &error),
            ELFF_XML_OK)
    << error.message;
  EXPECT_STREQ(value.data(), "flapping_fibre_output");
}

TEST_F(XMLCaseDocumentTest, TraversesOrderedElementChildrenOnly)
{
  elff_xml_error error = {};
  elff_xml_node_t root = nullptr, elff = nullptr, models = nullptr;
  elff_xml_node_t model = nullptr;
  size_t count = 0;
  const char* name = nullptr;

  ASSERT_EQ(elff_xml_document_root(document_, "Case", &root, &error),
            ELFF_XML_OK);
  ASSERT_EQ(elff_xml_node_unique_child(root, "ELFF", &elff, &error),
            ELFF_XML_OK);
  ASSERT_EQ(elff_xml_node_unique_child(elff, "Models", &models, &error),
            ELFF_XML_OK);
  ASSERT_EQ(elff_xml_node_element_child_count(models, &count, &error),
            ELFF_XML_OK);
  ASSERT_EQ(count, 1U);
  ASSERT_EQ(elff_xml_node_element_child_at(models, 0, &model, &error),
            ELFF_XML_OK);
  ASSERT_EQ(elff_xml_node_name(model, &name, &error), ELFF_XML_OK);
  EXPECT_STREQ(name, "EulerBeamHuang");
  EXPECT_EQ(elff_xml_node_element_child_at(models, 1, &model, &error),
            ELFF_XML_NOT_FOUND);
}

TEST_F(XMLCaseDocumentTest, LoadsAndValidatesOneOwnedSection)
{
  const std::string dtd(reinterpret_cast<const char*>(elff_embedded_elff_v2_dtd),
                        elff_embedded_elff_v2_dtd_size);
  elff_xml_document_t document = nullptr;
  elff_xml_node_t elff = nullptr;
  elff_xml_error error = {};
  const char* name = nullptr;

  ASSERT_EQ(elff_xml_case_load_section_file(ELFF_FLAPPING_FIBRE_CASE_PATH,
                                             "Case",
                                             "ELFF",
                                             dtd.data(),
                                             dtd.size(),
                                             &document,
                                             &elff,
                                             &error),
            ELFF_XML_OK)
    << error.message;
  ASSERT_NE(document, nullptr);
  ASSERT_NE(elff, nullptr);
  ASSERT_EQ(elff_xml_node_name(elff, &name, &error), ELFF_XML_OK);
  EXPECT_STREQ(name, "ELFF");

  elff_xml_document_destroy(document);
}

TEST(XMLCaseDocumentStandaloneTest, ReportsMalformedDocuments)
{
  constexpr char xml[] = "<Case><ELFF></Case>";
  elff_xml_document_t document = nullptr;
  elff_xml_error error = {};

  EXPECT_EQ(elff_xml_document_load_memory(
              xml, sizeof(xml) - 1, &document, &error),
            ELFF_XML_MALFORMED_DOCUMENT);
  EXPECT_EQ(document, nullptr);
  EXPECT_NE(error.message[0], '\0');
  EXPECT_GT(error.line, 0U);
}

TEST(XMLCaseDocumentStandaloneTest, RequiresUniqueDirectChildren)
{
  constexpr char xml[] = "<Case><ELFF/><ELFF/></Case>";
  elff_xml_document_t document = nullptr;
  elff_xml_node_t root = nullptr;
  elff_xml_node_t child = nullptr;
  elff_xml_error error = {};

  ASSERT_EQ(elff_xml_document_load_memory(
              xml, sizeof(xml) - 1, &document, &error),
            ELFF_XML_OK)
    << error.message;
  ASSERT_EQ(elff_xml_document_root(document, "Case", &root, &error),
            ELFF_XML_OK);
  EXPECT_EQ(elff_xml_node_unique_child(root, "ELFF", &child, &error),
            ELFF_XML_NOT_UNIQUE);
  EXPECT_EQ(child, nullptr);
  EXPECT_NE(error.message[0], '\0');

  elff_xml_document_destroy(document);
}

TEST(XMLCaseDocumentStandaloneTest, RejectsInvalidNumericAttributes)
{
  constexpr char xml[] = R"XML(
<Case>
  <Values integer="12x" real="nan"/>
</Case>
)XML";
  elff_xml_document_t document = nullptr;
  elff_xml_node_t root = nullptr;
  elff_xml_node_t values = nullptr;
  elff_xml_error error = {};
  int integer_value = 0;
  double real_value = 0.0;

  ASSERT_EQ(elff_xml_document_load_memory(
              xml, sizeof(xml) - 1, &document, &error),
            ELFF_XML_OK)
    << error.message;
  ASSERT_EQ(elff_xml_document_root(document, "Case", &root, &error),
            ELFF_XML_OK);
  ASSERT_EQ(elff_xml_node_unique_child(root, "Values", &values, &error),
            ELFF_XML_OK);

  EXPECT_EQ(elff_xml_node_attribute_int(
              values, "integer", &integer_value, &error),
            ELFF_XML_INVALID_VALUE);
  EXPECT_EQ(elff_xml_node_attribute_double(
              values, "real", &real_value, &error),
            ELFF_XML_INVALID_VALUE);

  elff_xml_document_destroy(document);
}

TEST(XMLCaseDocumentStandaloneTest, ReadsXYZWithOptionalZDefault)
{
  constexpr char xml[] = "<Case><Point x=\"1.5\" y=\"-2\"/></Case>";
  elff_xml_document_t document = nullptr;
  elff_xml_node_t root = nullptr;
  elff_xml_node_t point = nullptr;
  elff_xml_error error = {};
  double value[3] = {};

  ASSERT_EQ(elff_xml_document_load_memory(
              xml, sizeof(xml) - 1, &document, &error),
            ELFF_XML_OK)
    << error.message;
  ASSERT_EQ(elff_xml_document_root(document, "Case", &root, &error),
            ELFF_XML_OK);
  ASSERT_EQ(elff_xml_node_unique_child(root, "Point", &point, &error),
            ELFF_XML_OK);
  ASSERT_EQ(elff_xml_node_attribute_xyz(point, 3.25, value, &error),
            ELFF_XML_OK)
    << error.message;
  EXPECT_DOUBLE_EQ(value[0], 1.5);
  EXPECT_DOUBLE_EQ(value[1], -2.0);
  EXPECT_DOUBLE_EQ(value[2], 3.25);

  elff_xml_document_destroy(document);
}

TEST(XMLCaseDocumentStandaloneTest, DistinguishesMissingFiles)
{
  elff_xml_document_t document = nullptr;
  elff_xml_error error = {};

  EXPECT_EQ(elff_xml_document_load_file(
              "/definitely/not/an/elff/case.xml", &document, &error),
            ELFF_XML_IO_ERROR);
  EXPECT_EQ(document, nullptr);
  EXPECT_NE(error.message[0], '\0');
}

#else

TEST(XMLCaseDocumentTest, ReportsWhenLibXml2IsDisabled)
{
  constexpr char xml[] = "<Case/>";
  elff_xml_document_t document = nullptr;
  elff_xml_error error = {};

  EXPECT_EQ(elff_xml_document_load_memory(
              xml, sizeof(xml) - 1, &document, &error),
            ELFF_XML_UNAVAILABLE);
  EXPECT_EQ(document, nullptr);
  EXPECT_NE(error.message[0], '\0');
}

#endif

} // namespace
