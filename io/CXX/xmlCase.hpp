#pragma once

#include "elff/io/C/xmlCase.h"
#include <array>
#include <stdexcept>
#include <string>
#include <vector>

namespace ELFF::IO::CXX {

/** C++ diagnostic preserving the C XML status and source location. */
class XMLError : public std::runtime_error {
public:
  XMLError(elff_xml_status status, const elff_xml_error& diagnostic);
  elff_xml_status status() const noexcept { return status_; }
  const elff_xml_error& diagnostic() const noexcept { return diagnostic_; }
private:
  elff_xml_status status_;
  elff_xml_error diagnostic_;
};

/** Borrowed, read-only node. The source document must outlive every node view.
 * Operations delegate to the C API and throw XMLError on failure, including
 * ELFF_XML_UNAVAILABLE when XML support is disabled. No runtime/model policy. */
class XMLNode {
public:
  explicit XMLNode(elff_xml_node_t node) : node_(node) {}
  XMLNode child(const char* name) const;
  /** @brief Read finite numeric leaf text. */
  double text_number() const;
  /** @brief Return ordered element children with the requested name. */
  std::vector<XMLNode> children(const char* name) const;
  std::string name() const;
  size_t line() const;
  double number(const char* attribute) const;
  int integer(const char* attribute) const;
  std::string string(const char* attribute) const;
  std::array<double, 3> xyz(double default_z = 0.0) const;
  [[noreturn]] void invalid(const char* message) const;
private:
  elff_xml_node_t node_;
};

} // namespace ELFF::IO::CXX
