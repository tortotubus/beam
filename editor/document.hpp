#pragma once
#include "constraints.hpp"
#include <functional>
#include <optional>

namespace editor {

/**
 * @brief Own the XML tree and validate edits against the embedded DTD.
 */
class Document {
public:
  /**
   * @brief Load an existing XML file or start an unsaved document at a missing
   * path.
   */
  explicit Document(const std::string &filename, bool createIfMissing = true);

  /**
   * @brief Release the XML tree and its DTD.
   */
  ~Document();
  Document(const Document &) = delete;
  Document &operator=(const Document &) = delete;

  /**
   * @brief Borrow the tree for inspection; committed edits invalidate its
   * pointers.
   */
  xmlDocPtr xml() const { return tree; }

  /**
   * @brief Report whether an attribute still needs its required value filled.
   */
  bool isPending(xmlAttrPtr attribute) const {
    return pending.count(attribute) != 0;
  }

  /**
   * @brief Report whether the document has unsaved changes.
   */
  bool modified() const { return dirty; }
  std::string message;

  /**
   * @brief Check pending fields and full DTD validity, recording a status
   * message.
   */
  bool validate();

  /**
   * @brief Validate and write the document, creating a new file without
   * overwriting another.
   */
  bool save();

  /**
   * @brief List insertable elements at a child index that includes text and
   * comments.
   */
  std::vector<std::string> insertableElements(xmlNodePtr parent,
                                              std::size_t before) const;
  /**
   * @brief List constructible elements that may replace the selected element.
   */
  std::vector<std::string> replacementElements(xmlNodePtr node) const;

  /**
   * @brief Insert an element with required descendants and pending attribute
   * fields.
   */
  bool insertElement(xmlNodePtr parent, std::size_t before,
                     const std::string &name);

  /**
   * @brief Replace a non-root subtree with a newly constructed element.
   */
  bool replaceElement(xmlNodePtr node, const std::string &name);

  /**
   * @brief Delete a non-root subtree if the resulting document remains
   * permitted.
   */
  bool removeElement(xmlNodePtr node);

  /**
   * @brief Replace a text node with literal text, or remove it when the value
   * is empty.
   */
  bool setText(xmlNodePtr text, const std::string &value);

  /**
   * @brief Add literal text before a child index, or append it when no index
   * is supplied.
   */
  bool insertText(xmlNodePtr parent, const std::string &value,
                  std::optional<std::size_t> before = std::nullopt);

  /**
   * @brief Set a declared attribute value and resolve its pending field if
   * valid.
   */
  bool setAttribute(xmlNodePtr node, const std::string &name,
                    const std::string &value);

  /**
   * @brief Remove an explicit attribute unless it is required by the DTD.
   */
  bool removeAttribute(xmlNodePtr node, const std::string &name);

  /**
   * @brief Restore an attribute default or remove it; reject required and
   * fixed attributes.
   */
  bool resetAttribute(xmlNodePtr node, const std::string &name);

private:
  xmlDocPtr tree = nullptr;
  std::set<xmlAttrPtr> pending;
  bool dirty = false;
  std::string path;
  bool newFile = false;

  /**
   * @brief Apply an edit to a copy, committing or rolling back after
   * validation.
   */
  bool change(xmlNodePtr target, const std::function<bool(xmlNodePtr)> &action);
  /**
   * @brief Remove pending-attribute entries belonging to a subtree about to be
   * deleted.
   */
  void forget(xmlNodePtr node);
};
} // namespace editor
