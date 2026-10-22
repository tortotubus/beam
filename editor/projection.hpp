#pragma once
#include <libxml/tree.h>
#include <string>
#include <vector>

namespace editor {
class Document;

/**
 * @brief Escape control and non-ASCII bytes for predictable terminal display.
 */
std::string visible(const std::string &text);

enum class SyntaxColor {
  Delimiter = 1,
  TagName = 2,
  AttributeName = 3,
  Placeholder = 4
};

struct SyntaxSpan {
  std::size_t line;
  std::size_t column;
  std::size_t length;
  SyntaxColor color;
};

enum class Kind { Element, Attribute, AttributeGap, Text, TextGap, Gap };

struct Location {
  Kind kind;
  xmlNodePtr node;
  std::string attribute;
  // Physical child index, including whitespace and comments, for insertion.
  std::size_t before = 0;
  std::size_t line = 0;
  std::size_t column = 0;
  std::size_t length = 1;
};

/**
 * @brief Build display lines, cursor stops and colors; rebuild after committed
 * edits.
 */
struct Projection {
  /**
   * @brief Create an empty projection without accessing a document or
   * terminal.
   */
  Projection() = default;
  /**
   * @brief Render a document into display lines, borrowed node locations and
   * color spans.
   */
  explicit Projection(const Document &document);
  std::vector<std::string> lines;
  std::vector<Location> locations;
  std::vector<SyntaxSpan> styles;

private:
  /**
   * @brief Append an element subtree with its editable fields and permitted
   * insertion stops.
   */
  void renderNode(const Document &doc, xmlNodePtr node, int depth);
};
} // namespace editor
