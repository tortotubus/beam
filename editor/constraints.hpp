#pragma once
#include <libxml/tree.h>
#include <set>
#include <string>
#include <vector>

namespace editor {
/**
 * @brief Inspect DTD declarations and match ordered child sequences.
 */
class Constraints {
public:
  /**
   * @brief Borrow a document and its DTD declarations for constraint queries.
   */
  explicit Constraints(xmlDocPtr document) : doc(document) {}

  /**
   * @brief Find an element declaration, preferring the internal subset.
   */
  xmlElementPtr element(const xmlChar *name) const;

  /**
   * @brief List attribute declarations for an element, omitting shadowed
   * duplicates.
   */
  std::vector<xmlAttributePtr> attributes(const xmlChar *name) const;

  /**
   * @brief Return sorted, unique element names declared by the DTD.
   */
  std::vector<std::string> names() const;

  /**
   * @brief Match an ordered list of child-element names against the parent
   * content model.
   */
  bool accepts(xmlNodePtr parent,
               const std::vector<std::string> &children) const;
  /**
   * @brief List constructible elements permitted at an index in the child-
   * element sequence.
   */
  std::vector<std::string> insertions(xmlNodePtr parent,
                                      const std::vector<std::string> &children,
                                      std::size_t index) const;
  /**
   * @brief Report whether an element permits character data.
   */
  bool textAllowed(xmlNodePtr node) const;

  /**
   * @brief Report whether the DTD permits child elements, even if none exist.
   */
  bool childrenAllowed(xmlNodePtr node) const;

  /**
   * @brief Check child ordering and permitted text recursively, reporting the
   * first error.
   */
  bool structure(xmlNodePtr node, std::string &error) const;

  /**
   * @brief Build a detached minimal subtree and register its pending required
   * attributes.
   */
  xmlNodePtr create(const std::string &name,
                    std::set<xmlAttrPtr> &pending) const;

private:
  xmlDocPtr doc;
};

/**
 * @brief Find an explicit tree attribute without falling back to a DTD
 * default.
 */
xmlAttrPtr attributeNode(xmlNodePtr node, const xmlChar *name);

/**
 * @brief Copy a node or attribute content value into an owned string.
 */
std::string content(xmlNodePtr node);

/**
 * @brief Collect immediate child-element names in document order.
 */
std::vector<std::string> childNames(xmlNodePtr node);

} // namespace editor
