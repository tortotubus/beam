#include "document.hpp"
#include "schema.hpp"
#include <algorithm>
#include <cerrno>
#include <cstdarg>
#include <cstdio>
#include <cstring>
#include <libxml/chvalid.h>
#include <libxml/parser.h>
#include <libxml/valid.h>
#include <map>
#include <stdexcept>
#include <sys/stat.h>
#include <unistd.h>

namespace editor {
namespace {

struct ChildPosition {
  xmlNodePtr before;
  std::size_t elementIndex;
};
/**
 * @brief Resolve a physical child index and its corresponding element-only
 * index.
 */
std::optional<ChildPosition> childPosition(xmlNodePtr parent,
                                           std::size_t index) {
  if (!parent || parent->type != XML_ELEMENT_NODE)
    return std::nullopt;
  auto child = parent->children;
  std::size_t elementIndex = 0;
  while (index) {
    if (!child)
      return std::nullopt;
    if (child->type == XML_ELEMENT_NODE)
      ++elementIndex;
    child = child->next;
    --index;
  }
  return ChildPosition{child, elementIndex};
}

/**
 * @brief Count preceding siblings, including text and comment nodes.
 */
std::size_t siblingIndex(xmlNodePtr node) {
  std::size_t index = 0;
  for (auto previous = node->prev; previous; previous = previous->prev)
    ++index;
  return index;
}

/**
* @brief Check whether a node is an element beneath another element.
 */
bool replaceable(xmlNodePtr node) {
  return node && node->type == XML_ELEMENT_NODE && node->parent &&
         node->parent->type == XML_ELEMENT_NODE;
}

/**
 * @brief Find the DTD declaration for a named attribute of an element.
 */
xmlAttributePtr attributeDeclaration(xmlDocPtr tree, xmlNodePtr node,
                                     const std::string &name) {
  if (!node || node->type != XML_ELEMENT_NODE)
    return nullptr;
  for (auto declaration : Constraints(tree).attributes(node->name))
    if (name == reinterpret_cast<const char *>(declaration->name))
      return declaration;
  return nullptr;
}

/**
 * @brief Append a libxml2 diagnostic to the document status message.
 */
void diagnostic(void *context, const char *format, ...) {
  char buffer[2048];
  va_list args;
  va_start(args, format);
  vsnprintf(buffer, sizeof(buffer), format, args);
  va_end(args);
  auto &message = *static_cast<std::string *>(context);
  if (message.size() < 8192)
    message += buffer;
}

/**
 * @brief Parse the compiled schema and associate its declarations with a
 * document.
 */
xmlDtdPtr parseEmbeddedDtd(xmlDocPtr doc) {
  auto input = xmlParserInputBufferCreateMem(
      embeddedDtd, sizeof(embeddedDtd) - 1, XML_CHAR_ENCODING_UTF8);
  if (!input)
    throw std::runtime_error("Cannot allocate the embedded DTD input.");
  auto dtd = xmlIOParseDTD(nullptr, input, XML_CHAR_ENCODING_UTF8);
  if (!dtd)
    throw std::runtime_error("Cannot parse the embedded ELFF v2 DTD.");
  xmlSetTreeDoc(reinterpret_cast<xmlNodePtr>(dtd), doc);
  return dtd;
}

/**
 * @brief Require a Case root and validate the document against its attached
 * DTD.
 */
bool valid(xmlDocPtr doc, std::string &message) {
  xmlValidCtxt context{};
  context.userData = &message;
  context.error = diagnostic;
  context.warning = diagnostic;
  message.clear();
  if (!xmlDocGetRootElement(doc) ||
      !xmlStrEqual(xmlDocGetRootElement(doc)->name, BAD_CAST "Case")) {
    message = "The embedded ELFF v2 schema requires a <Case> root.";
    return false;
  }
  return xmlValidateDtd(&context, doc, doc->extSubset) == 1;
}

/**
 * @brief Reject namespaces and general entity references throughout a subtree.
 */
bool supported(xmlNodePtr node) {
  for (; node; node = node->next) {
    if (node->type == XML_ELEMENT_NODE) {
      if (node->ns || node->nsDef)
        return false;
      for (auto a = node->properties; a; a = a->next)
        if (a->ns || !supported(a->children))
          return false;
      if (!supported(node->children))
        return false;
    } else if (node->type == XML_ENTITY_REF_NODE)
      return false;
  }
  return true;
}

/**
 * @brief Check that a UTF-8 string contains only valid XML 1.0 characters.
 */
bool validText(const std::string &value) {
  for (std::size_t i = 0; i < value.size();) {
    int length = static_cast<int>(std::min<std::size_t>(4, value.size() - i));
    int code = xmlGetUTF8Char(
        reinterpret_cast<const unsigned char *>(value.data() + i), &length);
    if (length <= 0 || code < 0 || !xmlIsCharQ(code))
      return false;
    i += static_cast<std::size_t>(length);
  }
  return true;
}
} // namespace

Document::Document(const std::string &filename, bool createIfMissing) : path(filename) {
  struct stat info{};
  if (lstat(path.c_str(), &info) != 0) {
    if (errno != ENOENT || !createIfMissing)
      throw std::runtime_error("Cannot open " + path + ": " +
                               std::strerror(errno));
    auto slash = path.find_last_of('/');
    auto parent = slash == std::string::npos
                      ? "."
                      : (slash == 0 ? "/" : path.substr(0, slash));
    if (path.empty() || path.back() == '/' ||
        stat(parent.c_str(), &info) != 0 || !S_ISDIR(info.st_mode))
      throw std::runtime_error(
          "The destination's parent directory must exist: " + path);
    newFile = true;
    tree = xmlNewDoc(BAD_CAST "1.0");
  } else {
    if (!S_ISREG(info.st_mode))
      throw std::runtime_error("The document path must be a regular file: " +
                               path);
    // Do not load the document's external schema or materialize its defaults.
    tree = xmlReadFile(path.c_str(), nullptr, XML_PARSE_NONET);
  }
  if (!tree)
    throw std::runtime_error("Cannot parse or create XML: " + path);
  try {
    if (!newFile && !supported(xmlDocGetRootElement(tree)))
      throw std::runtime_error("Namespaces and general entity references are "
                               "not supported by this editor.");
    // The compiled schema is the only authority. Remove source declarations so
    // neither an internal subset nor a stale SYSTEM path can override it.
    if (tree->intSubset) {
      auto subset = tree->intSubset;
      xmlUnlinkNode(reinterpret_cast<xmlNodePtr>(subset));
      tree->intSubset = nullptr;
      xmlFreeDtd(subset);
    }
    if (tree->extSubset) {
      xmlFreeDtd(tree->extSubset);
      tree->extSubset = nullptr;
    }
    tree->extSubset = parseEmbeddedDtd(tree);
    if (newFile) {
      auto root = Constraints(tree).create("Case", pending);
      if (!root)
        throw std::runtime_error(
            "Cannot construct a minimal Case from the embedded DTD.");
      xmlDocSetRootElement(tree, root);
      dirty = true;
      message = "New ELFF v2 Case: fill required fields, then save to " + path;
    } else {
      // Apply defaults from the embedded schema, not from the source DOCTYPE.
      std::function<void(xmlNodePtr)> defaults = [&](auto node) {
        for (; node; node = node->next) {
          if (node->type != XML_ELEMENT_NODE)
            continue;
          for (auto declaration : Constraints(tree).attributes(node->name)) {
            if (!declaration->defaultValue)
              continue;
            bool present = false;
            for (auto attribute = node->properties; attribute;
                 attribute = attribute->next)
              if (xmlStrEqual(attribute->name, declaration->name))
                present = true;
            if (!present)
              xmlSetProp(node, declaration->name, declaration->defaultValue);
          }
          defaults(node->children);
        }
      };
      defaults(xmlDocGetRootElement(tree));
      if (!valid(tree, message))
        throw std::runtime_error(message);
      message = "Loaded using the embedded ELFF v2 DTD. ? for help.";
    }
  } catch (...) {
    xmlFreeDoc(tree);
    tree = nullptr;
    throw;
  }
}

Document::~Document() { xmlFreeDoc(tree); }

bool Document::validate() {
  if (!pending.empty()) {
    message = "Fill " + std::to_string(pending.size()) +
              " required attribute placeholder(s) before saving.";
    return false;
  }
  if (!valid(tree, message))
    return false;
  message = "Document is DTD-valid.";
  return true;
}

void Document::forget(xmlNodePtr node) {
  for (auto a = node->properties; a; a = a->next)
    pending.erase(a);
  for (auto c = node->children; c; c = c->next)
    if (c->type == XML_ELEMENT_NODE)
      forget(c);
}

bool Document::change(xmlNodePtr target,
                      const std::function<bool(xmlNodePtr)> &action) {
  auto original = tree;
  auto oldPending = pending;
  auto copy = xmlCopyDoc(tree, 1);
  if (!copy) {
    message = "Cannot copy document.";
    return false;
  }
  // Reparse the compiled schema for each copy. xmlCopyDoc omits the external
  // subset, and xmlCopyDtd does not preserve the element-to-attribute links
  // used to check required attributes in all supported libxml2 versions.
  try {
    if (copy->extSubset)
      xmlFreeDtd(copy->extSubset);
    copy->extSubset = nullptr;
    copy->extSubset = parseEmbeddedDtd(copy);
  } catch (...) {
    xmlFreeDoc(copy);
    throw;
  }
  std::map<xmlNodePtr, xmlNodePtr> nodes;
  std::set<xmlAttrPtr> copiedPending;
  std::function<void(xmlNodePtr, xmlNodePtr)> map = [&](auto a, auto b) {
    for (; a && b; a = a->next, b = b->next) {
      nodes[a] = b;
      if (a->type == XML_ELEMENT_NODE) {
        for (auto p = a->properties; p; p = p->next)
          if (pending.count(p))
            copiedPending.insert(attributeNode(b, p->name));
        map(a->children, b->children);
      }
    }
  };

  map(xmlDocGetRootElement(original), xmlDocGetRootElement(copy));
  tree = copy;
  pending = copiedPending;
  bool ok = false;

  try {
    ok = action(nodes.at(target));
    if (ok)
      ok = Constraints(tree).structure(xmlDocGetRootElement(tree), message);
    if (ok && pending.empty())
      ok = valid(tree, message);
  } catch (...) {
    tree = original;
    pending = oldPending;
    xmlFreeDoc(copy);
    throw;
  }
  if (!ok) {
    tree = original;
    pending = oldPending;
    xmlFreeDoc(copy);
    return false;
  }

  xmlFreeDoc(original);
  dirty = true;
  message =
      pending.empty()
          ? "Edit committed; document is valid."
          : "Incomplete: fill highlighted required attributes before saving.";
  return true;
}

std::vector<std::string>
Document::insertableElements(xmlNodePtr parent, std::size_t before) const {
  auto position = childPosition(parent, before);
  if (!position)
    return {};
  return Constraints(tree).insertions(parent, childNames(parent),
                                      position->elementIndex);
}

std::vector<std::string> Document::replacementElements(xmlNodePtr node) const {
  if (!replaceable(node))
    return {};
  auto position = childPosition(node->parent, siblingIndex(node));
  auto names = childNames(node->parent);
  names.erase(names.begin() + position->elementIndex);
  return Constraints(tree).insertions(node->parent, names,
                                      position->elementIndex);
}

bool Document::insertElement(xmlNodePtr parent, std::size_t before,
                             const std::string &name) {
  if (!childPosition(parent, before)) {
    message = "Invalid element insertion position.";
    return false;
  }
  return change(parent, [&](auto copy) {
    auto node = Constraints(tree).create(name, pending);
    if (!node) {
      message = "Required subtree exceeds the construction limit or has no "
                "finite construction.";
      return false;
    }
    auto position = childPosition(copy, before);
    if (position->before)
      xmlAddPrevSibling(position->before, node);
    else
      xmlAddChild(copy, node);
    return true;
  });
}

bool Document::replaceElement(xmlNodePtr node, const std::string &name) {
  if (!replaceable(node)) {
    message = "The document root cannot be removed or replaced.";
    return false;
  }
  return change(node, [&](auto copy) {
    auto replacement = Constraints(tree).create(name, pending);
    if (!replacement) {
      message = "Required subtree exceeds the construction limit or has no "
                "finite construction.";
      return false;
    }
    forget(copy);
    xmlReplaceNode(copy, replacement);
    xmlFreeNode(copy);
    return true;
  });
}

bool Document::removeElement(xmlNodePtr node) {
  if (!replaceable(node)) {
    message = "The document root cannot be removed or replaced.";
    return false;
  }
  auto position = childPosition(node->parent, siblingIndex(node));
  auto names = childNames(node->parent);
  names.erase(names.begin() + position->elementIndex);
  if (!Constraints(tree).accepts(node->parent, names)) {
    message = "Deletion would violate required children or cardinality.";
    return false;
  }
  return change(node, [&](auto copy) {
    forget(copy);
    xmlUnlinkNode(copy);
    xmlFreeNode(copy);
    return true;
  });
}

bool Document::setText(xmlNodePtr text, const std::string &value) {
  if (!text ||
      (text->type != XML_TEXT_NODE && text->type != XML_CDATA_SECTION_NODE) ||
      !text->parent || !Constraints(tree).textAllowed(text->parent)) {
    message = "Select an editable text node.";
    return false;
  }
  if (!validText(value)) {
    message = "Text must contain valid XML 1.0 characters encoded as UTF-8.";
    return false;
  }
  return change(text, [&](auto copy) {
    // Setting content on a text node treats the input as literal text.
    if (value.empty()) {
      xmlUnlinkNode(copy);
      xmlFreeNode(copy);
    } else
      xmlNodeSetContent(copy, BAD_CAST value.c_str());
    return true;
  });
}

bool Document::insertText(xmlNodePtr parent, const std::string &value,
                          std::optional<std::size_t> before) {
  if (!parent || parent->type != XML_ELEMENT_NODE ||
      !Constraints(tree).textAllowed(parent)) {
    message = "This element does not allow character data.";
    return false;
  }
  if (before && !childPosition(parent, *before)) {
    message = "Invalid text insertion position.";
    return false;
  }
  if (value.empty() && !parent->children)
    return true;
  if (!validText(value)) {
    message = "Text must contain valid XML 1.0 characters encoded as UTF-8.";
    return false;
  }
  return change(parent, [&](auto copy) {
    auto text = xmlNewDocText(tree, BAD_CAST value.c_str());
    auto next = before ? childPosition(copy, *before)->before : nullptr;
    if (next)
      xmlAddPrevSibling(next, text);
    else
      xmlAddChild(copy, text);
    return true;
  });
}

bool Document::setAttribute(xmlNodePtr node, const std::string &name,
                            const std::string &value) {
  if (!attributeDeclaration(tree, node, name)) {
    message = "Undeclared attribute.";
    return false;
  }
  if (!validText(value)) {
    message = "Value must contain valid XML 1.0 characters encoded as UTF-8.";
    return false;
  }
  return change(node, [&](auto copy) {
    auto old = attributeNode(copy, BAD_CAST name.c_str());
    pending.erase(old);
    auto attr = xmlSetProp(copy, BAD_CAST name.c_str(), BAD_CAST value.c_str());
    xmlValidCtxt context{};
    context.userData = &message;
    context.error = diagnostic;
    context.warning = diagnostic;
    message.clear();
    return xmlValidateOneAttribute(&context, tree, copy, attr,
                                   BAD_CAST value.c_str()) == 1;
  });
}

bool Document::removeAttribute(xmlNodePtr node, const std::string &name) {
  auto declaration = attributeDeclaration(tree, node, name);
  if (!declaration) {
    message = "Attribute has no declaration.";
    return false;
  }
  if (declaration->def == XML_ATTRIBUTE_REQUIRED) {
    message = "Required attributes cannot be removed.";
    return false;
  }
  return change(node, [&](auto copy) {
    auto attr = attributeNode(copy, BAD_CAST name.c_str());
    if (attr) {
      pending.erase(attr);
      xmlRemoveProp(attr);
    }
    return true;
  });
}

bool Document::resetAttribute(xmlNodePtr node, const std::string &name) {
  auto declaration = attributeDeclaration(tree, node, name);
  if (!declaration) {
    message = "Attribute has no declaration.";
    return false;
  }
  if (declaration->def == XML_ATTRIBUTE_REQUIRED) {
    message = "Required attributes cannot be removed.";
    return false;
  }
  if (declaration->def == XML_ATTRIBUTE_FIXED) {
    message = "This attribute is fixed by the DTD.";
    return false;
  }
  if (declaration->defaultValue)
    return setAttribute(
        node, name, reinterpret_cast<const char *>(declaration->defaultValue));
  return removeAttribute(node, name);
}

bool Document::save() {
  if (!validate())
    return false;
  struct stat info{};
  if (newFile) {
    if (lstat(path.c_str(), &info) == 0 || errno != ENOENT) {
      message = "The new document destination is no longer available; it will "
                "not be overwritten.";
      return false;
    }
    info.st_mode = 0600;
  } else if (lstat(path.c_str(), &info) != 0 || !S_ISREG(info.st_mode)) {
    message = "Save requires a regular file (symlinks are not replaced).";
    return false;
  }
  std::string pattern = path + ".tmp.XXXXXX";
  std::vector<char> temp(pattern.begin(), pattern.end());
  temp.push_back('\0');
  int fd = mkstemp(temp.data());
  if (fd < 0) {
    message = std::strerror(errno);
    return false;
  }
  xmlChar *data = nullptr;
  int size = 0;
  // Do not pretty-print on disk: indentation would change mixed-content text.
  xmlDocDumpMemoryEnc(tree, &data, &size, "UTF-8");
  bool ok = data && size > 0 && fchmod(fd, info.st_mode & 0777) == 0;
  int written = 0;
  while (ok && written < size) {
    auto count =
        write(fd, data + written, static_cast<std::size_t>(size - written));
    if (count < 0 && errno == EINTR)
      continue;
    if (count <= 0) {
      ok = false;
      break;
    }
    written += static_cast<int>(count);
  }
  xmlFree(data);
  if (ok && fsync(fd) != 0)
    ok = false;
  if (close(fd) != 0)
    ok = false;
  if (ok) {
    // Linking the completed temporary file publishes a new path atomically,
    // refusing to overwrite a file created since this editing session began.
    if (newFile) {
      if (link(temp.data(), path.c_str()) != 0)
        ok = false;
      else
        unlink(temp.data());
    } else if (rename(temp.data(), path.c_str()) != 0)
      ok = false;
  }
  if (!ok) {
    int error = errno;
    unlink(temp.data());
    message = "Save failed: " + std::string(std::strerror(error));
    return false;
  }
  newFile = false;
  dirty = false;
  message = "Saved " + path;
  return true;
}

} // namespace editor
