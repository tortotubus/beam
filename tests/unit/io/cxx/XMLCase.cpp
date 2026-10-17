#include <gtest/gtest.h>
#include "elff/config/_config.hpp"
#include "elff/io/CXX/xmlCase.hpp"
#include <cstring>
#include <memory>

using namespace ELFF::IO::CXX;

#ifdef ELFF_USE_LIBXML2
TEST(XMLNodeTest, ReadsTypedValuesAndOwnsReturnedStrings)
{
  const char* xml = "<Root><Point x='1.5' y='-2' count='3' label='hello'/></Root>";
  elff_xml_document_t raw = nullptr;
  ASSERT_EQ(elff_xml_document_load_memory(xml, std::strlen(xml), &raw, nullptr), ELFF_XML_OK);
  std::unique_ptr<elff_xml_document_handle, decltype(&elff_xml_document_destroy)>
    document(raw, elff_xml_document_destroy);
  elff_xml_node_t root = nullptr;
  ASSERT_EQ(elff_xml_document_root(raw, "Root", &root, nullptr), ELFF_XML_OK);
  auto point = XMLNode(root).child("Point");
  EXPECT_EQ(point.name(), "Point");
  EXPECT_EQ(point.integer("count"), 3);
  EXPECT_DOUBLE_EQ(point.number("x"), 1.5);
  EXPECT_EQ(point.xyz(4), (std::array<double, 3>{1.5, -2, 4}));
  const auto label = point.string("label");
  document.reset();
  EXPECT_EQ(label, "hello");
}

TEST(XMLNodeTest, PreservesStatusesAndSourceLines)
{
  const char* xml = "<Root>\n<Point count='bad'/>\n<Point/>\n</Root>";
  elff_xml_document_t raw = nullptr;
  ASSERT_EQ(elff_xml_document_load_memory(xml, std::strlen(xml), &raw, nullptr), ELFF_XML_OK);
  std::unique_ptr<elff_xml_document_handle, decltype(&elff_xml_document_destroy)>
    document(raw, elff_xml_document_destroy);
  elff_xml_node_t root = nullptr, point = nullptr;
  ASSERT_EQ(elff_xml_document_root(raw, "Root", &root, nullptr), ELFF_XML_OK);
  ASSERT_EQ(elff_xml_node_element_child_at(root, 0, &point, nullptr), ELFF_XML_OK);
  try { XMLNode(root).child("Point"); FAIL(); }
  catch (const XMLError& error) { EXPECT_EQ(error.status(), ELFF_XML_NOT_UNIQUE); }
  try { XMLNode(root).child("Missing"); FAIL(); }
  catch (const XMLError& error) { EXPECT_EQ(error.status(), ELFF_XML_NOT_FOUND); }
  try { XMLNode(point).integer("count"); FAIL(); }
  catch (const XMLError& error) {
    EXPECT_EQ(error.status(), ELFF_XML_INVALID_VALUE);
    EXPECT_EQ(error.diagnostic().line, 2u);
  }
  try { XMLNode(point).string("missing"); FAIL(); }
  catch (const XMLError& error) { EXPECT_EQ(error.status(), ELFF_XML_NOT_FOUND); }
  try { XMLNode(point).invalid("Model constraint failed"); FAIL(); }
  catch (const XMLError& error) {
    EXPECT_EQ(error.status(), ELFF_XML_INVALID_VALUE);
    EXPECT_EQ(error.diagnostic().line, 2u);
    EXPECT_STREQ(error.what(), "Model constraint failed");
  }
}
#else
TEST(XMLNodeTest, PreservesUnavailableStatus)
{
  // Disabled C stubs require a non-null handle but never dereference it.
  const char sentinel = 0;
  const auto node = reinterpret_cast<elff_xml_node_t>(&sentinel);
  try { XMLNode(node).number("length"); FAIL(); }
  catch (const XMLError& error) { EXPECT_EQ(error.status(), ELFF_XML_UNAVAILABLE); }
}
#endif
