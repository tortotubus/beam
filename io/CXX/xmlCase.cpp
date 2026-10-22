#include "elff/io/CXX/xmlCase.hpp"
#include <vector>

namespace ELFF::IO::CXX {
namespace {
void check(elff_xml_status status, const elff_xml_error& diagnostic)
{
  if (status != ELFF_XML_OK)
    throw XMLError(status, diagnostic);
}
}

XMLError::XMLError(elff_xml_status status, const elff_xml_error& diagnostic)
  : std::runtime_error(diagnostic.message[0] ? diagnostic.message
                                            : elff_xml_status_string(status)),
    status_(status), diagnostic_(diagnostic)
{}

XMLNode XMLNode::child(const char* name) const
{
  elff_xml_error error{};
  elff_xml_node_t result = nullptr;
  check(elff_xml_node_unique_child(node_, name, &result, &error), error);
  return XMLNode(result);
}

double XMLNode::text_number() const
{
  elff_xml_error error{};
  double result = 0;
  check(elff_xml_node_text_double(node_, &result, &error), error);
  return result;
}

std::vector<XMLNode> XMLNode::children(const char* name) const
{
  elff_xml_error error{};
  size_t count = 0;
  check(elff_xml_node_element_child_count(node_, &count, &error), error);
  std::vector<XMLNode> result;
  for (size_t i = 0; i < count; ++i) {
    elff_xml_node_t child = nullptr;
    check(elff_xml_node_element_child_at(node_, i, &child, &error), error);
    XMLNode node(child);
    if (node.name() == name) result.push_back(node);
  }
  return result;
}

std::string XMLNode::name() const
{
  elff_xml_error error{};
  const char* result = nullptr;
  check(elff_xml_node_name(node_, &result, &error), error);
  return result;
}

size_t XMLNode::line() const
{
  elff_xml_error error{};
  size_t result = 0;
  check(elff_xml_node_line(node_, &result, &error), error);
  return result;
}

double XMLNode::number(const char* attribute) const
{
  elff_xml_error error{};
  double result = 0;
  check(elff_xml_node_attribute_double(node_, attribute, &result, &error), error);
  return result;
}

int XMLNode::integer(const char* attribute) const
{
  elff_xml_error error{};
  int result = 0;
  check(elff_xml_node_attribute_int(node_, attribute, &result, &error), error);
  return result;
}

std::string XMLNode::string(const char* attribute) const
{
  elff_xml_error error{};
  size_t capacity = 0;
  const auto status = elff_xml_node_attribute_string(
    node_, attribute, nullptr, 0, &capacity, &error);
  if (status != ELFF_XML_BUFFER_TOO_SMALL)
    throw XMLError(status, error);
  std::vector<char> buffer(capacity);
  check(elff_xml_node_attribute_string(node_, attribute, buffer.data(),
                                      buffer.size(), nullptr, &error), error);
  return buffer.data();
}

std::array<double, 3> XMLNode::xyz(double default_z) const
{
  elff_xml_error error{};
  std::array<double, 3> result{};
  check(elff_xml_node_attribute_xyz(node_, default_z, result.data(), &error), error);
  return result;
}

void XMLNode::invalid(const char* message) const
{
  elff_xml_error error{};
  elff_xml_error_set(&error, line(), 0, message);
  throw XMLError(ELFF_XML_INVALID_VALUE, error);
}
} // namespace ELFF::IO::CXX
