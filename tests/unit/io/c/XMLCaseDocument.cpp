#include <gtest/gtest.h>

#include <elff/config/_config.hpp>
#include <elff/io/C/xmlCaseDocument.h>
#include <elff/io/C/xmlDTDValidator.h>

#include <fstream>
#include <iterator>
#include <string>
#include <vector>

namespace {

std::string read_text_file(const char* path)
{
  std::ifstream input(path, std::ios::binary);
  return std::string(std::istreambuf_iterator<char>(input),
                     std::istreambuf_iterator<char>());
}

#ifdef ELFF_USE_LIBXML2

class XMLCaseDocumentTest : public ::testing::Test
{
protected:
  void SetUp() override
  {
    elff_xml_error error = {};
    ASSERT_EQ(elff_xml_document_load_file(
                ELFF_FLAPPING_FIBRE_CASE_PATH, &document_, &error),
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
  EXPECT_EQ(format_version, 1);

  ASSERT_EQ(elff_xml_node_unique_child(
              root, "Basilisk", &basilisk, &error),
            ELFF_XML_OK)
    << error.message;
  ASSERT_EQ(elff_xml_node_unique_child(
              basilisk, "Domain", &domain, &error),
            ELFF_XML_OK)
    << error.message;
  ASSERT_EQ(elff_xml_node_unique_child(basilisk, "Grid", &grid, &error),
            ELFF_XML_OK)
    << error.message;

  EXPECT_EQ(elff_xml_node_attribute_double(
              domain, "size", &domain_size, &error),
            ELFF_XML_OK)
    << error.message;
  EXPECT_DOUBLE_EQ(domain_size, 16.0);
  EXPECT_EQ(elff_xml_node_attribute_int(
              grid, "minLevel", &min_level, &error),
            ELFF_XML_OK)
    << error.message;
  EXPECT_EQ(min_level, 6);
}

TEST_F(XMLCaseDocumentTest, ValidatesOwnedSectionsIndependently)
{
  const std::string elff_dtd = read_text_file(ELFF_CASE_DTD_PATH);
  const std::string basilisk_dtd = read_text_file(ELFF_BASILISK_DTD_PATH);
  elff_xml_error error = {};
  elff_xml_node_t root = nullptr;
  elff_xml_node_t elff = nullptr;
  elff_xml_node_t basilisk = nullptr;
  char diagnostic[512] = {};

  ASSERT_FALSE(elff_dtd.empty());
  ASSERT_FALSE(basilisk_dtd.empty());
  ASSERT_EQ(elff_xml_document_root(document_, "Case", &root, &error),
            ELFF_XML_OK)
    << error.message;
  ASSERT_EQ(elff_xml_node_unique_child(root, "ELFF", &elff, &error),
            ELFF_XML_OK)
    << error.message;
  ASSERT_EQ(elff_xml_node_unique_child(
              root, "Basilisk", &basilisk, &error),
            ELFF_XML_OK)
    << error.message;

  EXPECT_EQ(elff_xml_dtd_validate_node(elff,
                                        elff_dtd.data(),
                                        elff_dtd.size(),
                                        diagnostic,
                                        sizeof(diagnostic)),
            ELFF_XML_DTD_VALID)
    << diagnostic;
  EXPECT_EQ(elff_xml_dtd_validate_node(basilisk,
                                        basilisk_dtd.data(),
                                        basilisk_dtd.size(),
                                        diagnostic,
                                        sizeof(diagnostic)),
            ELFF_XML_DTD_VALID)
    << diagnostic;
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
