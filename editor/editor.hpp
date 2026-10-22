#pragma once
#include "cursor.hpp"
#include "document.hpp"
#include "terminal.hpp"

namespace editor {
/** @brief Coordinate terminal input, projection navigation and document
 * operations. */
class Editor {
public:
  /**
   * @brief Borrow a document and initialize the terminal interface.
   */
  explicit Editor(Document &document) : doc(document) {}
  /**
   * @brief Display the document and dispatch navigation and editing commands
   * until quit.
   */
  void run();

private:
  Document &doc;
  Terminal terminal;
  Cursor cursor;
  /**
   * @brief Handle element insertion, replacement or deletion at the selected
   * location.
   */
  void structural(const Location &location, int command);
  /**
   * @brief Edit a named attribute, or choose a missing attribute to add.
   */
  void attribute(xmlNodePtr node, const std::string &name = "");
  /**
   * @brief Prompt for text and apply it to the selected node or empty text
   * field.
   */
  void text(const Location &location);
};
} // namespace editor
