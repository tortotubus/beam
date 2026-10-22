#include "projection.hpp"
#include "document.hpp"
#include <algorithm>
#include <cstdio>
#include <libxml/parser.h>

namespace editor {
std::string visible(const std::string &text) {
  std::string out;
  for (unsigned char c : text) {
    if (c == '\n')
      out += "\\n";
    else if (c == '\r')
      out += "\\r";
    else if (c == '\t')
      out += "\\t";
    else if (c < 32 || c >= 127) {
      char buffer[5];
      snprintf(buffer, sizeof(buffer), "\\x%02x", c);
      out += buffer;
    } else
      out += static_cast<char>(c);
  }
  return out;
}
namespace {
/**
 * @brief Escape XML-sensitive characters and make the result safe for terminal
 * display.
 */
std::string escaped(xmlDocPtr doc, const std::string &value) {
  auto result = xmlEncodeSpecialChars(doc, BAD_CAST value.c_str());
  std::string text = result ? reinterpret_cast<const char *>(result) : "";
  xmlFree(result);
  return visible(text);
}
} // namespace

void Projection::renderNode(const Document &doc, xmlNodePtr node, int depth) {
  auto pad =
      std::string(static_cast<std::size_t>(std::min(depth, 30)) * 2, ' ');
  std::string name(reinterpret_cast<const char *>(node->name));
  const auto line = lines.size();
  locations.push_back(
      {Kind::Element, node, "", 0, line, pad.size() + 1, name.size()});
  styles.push_back({line, pad.size(), 1, SyntaxColor::Delimiter});
  styles.push_back({line, pad.size() + 1, name.size(), SyntaxColor::TagName});
  std::string opening = pad + "<" + name;

  for (auto a = node->properties; a; a = a->next) {
    std::string attributeName(reinterpret_cast<const char *>(a->name));
    styles.push_back({line, opening.size() + 1, attributeName.size(),
                      SyntaxColor::AttributeName});
    opening += " " + attributeName + "=\"";
    auto value =
        doc.isPending(a)
            ? " "
            : escaped(doc.xml(), content(reinterpret_cast<xmlNodePtr>(a)));
    // An empty value selects its closing quote so it still has a visible
    // cursor.
    locations.push_back({Kind::Attribute, node, attributeName, 0, line,
                         opening.size(),
                         std::max<std::size_t>(1, value.size())});
    if (doc.isPending(a))
      styles.push_back(
          {line, opening.size(), value.size(), SyntaxColor::Placeholder});
    opening += value + "\"";
  }
  Constraints constraints(doc.xml());
  auto attributes = constraints.attributes(node->name);
  bool missingAttributes =
      std::any_of(attributes.begin(), attributes.end(), [&](auto declaration) {
        return !attributeNode(node, declaration->name);
      });
  if (missingAttributes) {
    locations.push_back(
        {Kind::AttributeGap, node, "", 0, line, opening.size(), 1});
    opening += " ";
  }
  auto names = childNames(node);
  std::size_t elementIndex = 0;
  // A gap shares the trailing blank cell of the preceding rendered line.
  // The same query drives both cursor stops and the insertion menu.
  auto gap = [&](std::size_t index) {
    if (constraints.insertions(node, names, elementIndex).empty())
      return;
    auto &preceding = lines.back();
    locations.push_back(
        {Kind::Gap, node, "", index, lines.size() - 1, preceding.size(), 1});
    preceding += " ";
  };
  if (!node->children) {
    if (constraints.childrenAllowed(node) || constraints.textAllowed(node)) {
      styles.push_back({line, opening.size(), 1, SyntaxColor::Delimiter});
      lines.push_back(opening + ">");
      // Put the insertion stop inside the paired tags, before the closing tag.
      gap(0);
      auto &shown = lines.back();
      if (constraints.textAllowed(node)) {
        locations.push_back(
            {Kind::TextGap, node, "", 0, line, shown.size(), 1});
        styles.push_back({line, shown.size(), 1, SyntaxColor::Placeholder});
        shown += " ";
      }
      styles.push_back({line, shown.size(), 2, SyntaxColor::Delimiter});
      styles.push_back(
          {line, shown.size() + 2, name.size(), SyntaxColor::TagName});
      styles.push_back(
          {line, shown.size() + 2 + name.size(), 1, SyntaxColor::Delimiter});
      shown += "</" + name + ">";
    } else {
      if (!missingAttributes)
        opening += " ";
      lines.push_back(opening + "/>");
      styles.push_back(
          {line, lines.back().size() - 2, 2, SyntaxColor::Delimiter});
    }
    return;
  }
  bool textOnly = constraints.textAllowed(node) && names.empty() &&
                  constraints.insertions(node, names, 0).empty();
  for (auto child = node->children; child; child = child->next)
    if (child->type != XML_TEXT_NODE && child->type != XML_CDATA_SECTION_NODE)
      textOnly = false;
  styles.push_back({line, opening.size(), 1, SyntaxColor::Delimiter});
  if (textOnly) {
    opening += ">";
    for (auto child = node->children; child; child = child->next) {
      auto value = escaped(doc.xml(), content(child));
      locations.push_back({Kind::Text, child, "", 0, line, opening.size(),
                           std::max<std::size_t>(1, value.size())});
      if (value.empty())
        styles.push_back({line, opening.size(), 1, SyntaxColor::Placeholder});
      opening += value.empty() ? " " : value;
    }
    styles.push_back({line, opening.size(), 2, SyntaxColor::Delimiter});
    styles.push_back(
        {line, opening.size() + 2, name.size(), SyntaxColor::TagName});
    styles.push_back(
        {line, opening.size() + 2 + name.size(), 1, SyntaxColor::Delimiter});
    lines.push_back(opening + "</" + name + ">");
    return;
  }
  lines.push_back(opening + ">");
  bool hasText = false;
  for (auto child = node->children; child; child = child->next)
    hasText |=
        child->type == XML_TEXT_NODE || child->type == XML_CDATA_SECTION_NODE;
  if (constraints.textAllowed(node) && !hasText) {
    locations.push_back(
        {Kind::TextGap, node, "", 0, lines.size(), pad.size() + 2, 1});
    styles.push_back(
        {lines.size(), pad.size() + 2, 1, SyntaxColor::Placeholder});
    lines.push_back(pad + "   ");
  }
  std::size_t index = 0;
  for (auto child = node->children; child; child = child->next, ++index) {
    bool isText =
        child->type == XML_TEXT_NODE || child->type == XML_CDATA_SECTION_NODE;
    // Formatting whitespace in element-only content is not an editing stop.
    // Keep its actual sibling index so insertion still targets the XML tree.
    if (isText && !constraints.textAllowed(node))
      continue;
    gap(index);
    if (child->type == XML_ELEMENT_NODE) {
      renderNode(doc, child, depth + 1);
      ++elementIndex;
    } else if (isText) {
      auto value = escaped(doc.xml(), content(child));
      locations.push_back({Kind::Text, child, "", 0, lines.size(),
                           pad.size() + 2,
                           std::max<std::size_t>(1, value.size())});
      if (value.empty())
        styles.push_back(
            {lines.size(), pad.size() + 2, 1, SyntaxColor::Placeholder});
      lines.push_back(pad + "  " + (value.empty() ? " " : value));
    } else {
      auto buffer = xmlBufferCreate();
      xmlNodeDump(buffer, doc.xml(), child, 0, 0);
      lines.push_back(
          pad + "  " +
          visible(reinterpret_cast<const char *>(xmlBufferContent(buffer))));
      xmlBufferFree(buffer);
    }
  }
  gap(index);
  styles.push_back({lines.size(), pad.size(), 2, SyntaxColor::Delimiter});
  styles.push_back(
      {lines.size(), pad.size() + 2, name.size(), SyntaxColor::TagName});
  styles.push_back(
      {lines.size(), pad.size() + 2 + name.size(), 1, SyntaxColor::Delimiter});
  lines.push_back(pad + "</" + name + ">");
}

Projection::Projection(const Document &doc) {
  lines.push_back("<?xml version=\"1.0\"?>");
  for (auto n = doc.xml()->children; n; n = n->next) {
    if (n->type == XML_ELEMENT_NODE)
      renderNode(doc, n, 0);
    else {
      auto buffer = xmlBufferCreate();
      xmlNodeDump(buffer, doc.xml(), n, 0, 0);
      lines.push_back(
          visible(reinterpret_cast<const char *>(xmlBufferContent(buffer))));
      xmlBufferFree(buffer);
    }
  }
}
} // namespace editor
