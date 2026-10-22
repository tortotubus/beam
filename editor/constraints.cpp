#include "constraints.hpp"
#include <algorithm>
#include <functional>
#include <libxml/hash.h>
#include <libxml/valid.h>
#include <map>
#include <tuple>

namespace editor {

namespace {
using Positions = std::set<std::size_t>;

class Matcher {
  const std::vector<std::string> &input;
  std::map<std::tuple<xmlElementContentPtr, std::size_t, bool>, Positions>
      cache;

public:
  /**
   * @brief Borrow the child-name sequence to be matched.
   */
  explicit Matcher(const std::vector<std::string> &names) : input(names) {}

  /**
   * @brief Return reachable end positions for a model, optionally ignoring its
   * occurrence modifier.
   */
  Positions match(xmlElementContentPtr c, std::size_t start,
                  bool once = false) {
    if (!c)
      return {start};
    auto key = std::make_tuple(c, start, once);
    auto found = cache.find(key);
    if (found != cache.end())
      return found->second;
    Positions result;
    if (!once && c->ocur != XML_ELEMENT_CONTENT_ONCE) {
      result = match(c, start, true);
      if (c->ocur == XML_ELEMENT_CONTENT_OPT ||
          c->ocur == XML_ELEMENT_CONTENT_MULT)
        result.insert(start);
      if (c->ocur == XML_ELEMENT_CONTENT_MULT ||
          c->ocur == XML_ELEMENT_CONTENT_PLUS) {
        std::vector<std::size_t> work(result.begin(), result.end());
        for (std::size_t i = 0; i < work.size(); ++i)
          for (auto end : match(c, work[i], true))
            if (result.insert(end).second)
              work.push_back(end);
      }
    } else {
      switch (c->type) {
      case XML_ELEMENT_CONTENT_PCDATA:
        result.insert(start);
        break;
      case XML_ELEMENT_CONTENT_ELEMENT:
        if (start < input.size() &&
            input[start] == reinterpret_cast<const char *>(c->name))
          result.insert(start + 1);
        break;
      case XML_ELEMENT_CONTENT_OR:
        result = match(c->c1, start);
        {
          auto right = match(c->c2, start);
          result.insert(right.begin(), right.end());
        }
        break;
      case XML_ELEMENT_CONTENT_SEQ:
        for (auto middle : match(c->c1, start)) {
          auto right = match(c->c2, middle);
          result.insert(right.begin(), right.end());
        }
        break;
      }
    }
    cache[key] = result;
    return result;
  }
};
} // namespace

xmlAttrPtr attributeNode(xmlNodePtr node, const xmlChar *name) {
  for (auto attribute = node->properties; attribute;
       attribute = attribute->next)
    if (xmlStrEqual(attribute->name, name))
      return attribute;
  return nullptr;
}

std::string content(xmlNodePtr node) {
  xmlChar *value = xmlNodeGetContent(node);
  std::string result = value ? reinterpret_cast<const char *>(value) : "";
  xmlFree(value);
  return result;
}

std::vector<std::string> childNames(xmlNodePtr node) {
  std::vector<std::string> result;
  for (auto c = node->children; c; c = c->next)
    if (c->type == XML_ELEMENT_NODE)
      result.emplace_back(reinterpret_cast<const char *>(c->name));
  return result;
}

xmlElementPtr Constraints::element(const xmlChar *name) const {
  auto e =
      doc->intSubset ? xmlGetDtdElementDesc(doc->intSubset, name) : nullptr;
  return e ? e
           : (doc->extSubset ? xmlGetDtdElementDesc(doc->extSubset, name)
                             : nullptr);
}

std::vector<std::string> Constraints::names() const {
  std::set<std::string> unique;
  for (auto dtd : {doc->intSubset, doc->extSubset})
    if (dtd && dtd->elements)
      xmlHashScan(
          static_cast<xmlHashTablePtr>(dtd->elements),
          [](void *, void *data, const xmlChar *name) {
            static_cast<std::set<std::string> *>(data)->insert(
                reinterpret_cast<const char *>(name));
          },
          &unique);
  return {unique.begin(), unique.end()};
}

std::vector<xmlAttributePtr>
Constraints::attributes(const xmlChar *name) const {
  std::vector<xmlAttributePtr> result;
  for (auto dtd : {doc->intSubset, doc->extSubset}) {
    if (!dtd)
      continue;
    for (auto n = dtd->children; n; n = n->next) {
      if (n->type != XML_ATTRIBUTE_DECL)
        continue;
      auto a = reinterpret_cast<xmlAttributePtr>(n);
      if (!xmlStrEqual(a->elem, name))
        continue;
      if (std::none_of(result.begin(), result.end(),
                       [a](auto b) { return xmlStrEqual(a->name, b->name); }))
        result.push_back(a);
    }
  }
  return result;
}

bool Constraints::accepts(xmlNodePtr parent,
                          const std::vector<std::string> &children) const {
  auto e = element(parent->name);
  if (!e)
    return false;
  if (e->etype == XML_ELEMENT_TYPE_ANY)
    return true;
  if (e->etype == XML_ELEMENT_TYPE_EMPTY)
    return children.empty();
  return Matcher(children).match(e->content, 0).count(children.size()) != 0;
}

std::vector<std::string>
Constraints::insertions(xmlNodePtr parent,
                        const std::vector<std::string> &children,
                        std::size_t index) const {
  std::vector<std::string> result;
  for (const auto &name : names()) {
    auto candidate = children;
    candidate.insert(candidate.begin() + index, name);
    if (!accepts(parent, candidate))
      continue;
    std::set<xmlAttrPtr> temporary;
    auto node = create(name, temporary);
    if (node) {
      result.push_back(name);
      xmlFreeNode(node);
    }
  }
  return result;
}

bool Constraints::textAllowed(xmlNodePtr node) const {
  auto e = element(node->name);
  return e && (e->etype == XML_ELEMENT_TYPE_MIXED ||
               e->etype == XML_ELEMENT_TYPE_ANY);
}

bool Constraints::childrenAllowed(xmlNodePtr node) const {
  auto declaration = element(node->name);
  if (!declaration)
    return false;
  if (declaration->etype == XML_ELEMENT_TYPE_ANY)
    return true;
  std::vector<xmlElementContentPtr> pending{declaration->content};
  while (!pending.empty()) {
    auto model = pending.back();
    pending.pop_back();
    if (!model)
      continue;
    if (model->type == XML_ELEMENT_CONTENT_ELEMENT)
      return true;
    pending.push_back(model->c1);
    pending.push_back(model->c2);
  }
  return false;
}

bool Constraints::structure(xmlNodePtr node, std::string &error) const {
  if (!accepts(node, childNames(node))) {
    error = "Children of <" +
            std::string(reinterpret_cast<const char *>(node->name)) +
            "> violate its DTD content model.";
    return false;
  }
  auto e = element(node->name);
  for (auto c = node->children; c; c = c->next) {
    if (c->type == XML_ELEMENT_NODE && !structure(c, error))
      return false;
    if ((c->type == XML_TEXT_NODE || c->type == XML_CDATA_SECTION_NODE) &&
        !textAllowed(node) &&
        (e->etype == XML_ELEMENT_TYPE_EMPTY ||
         content(c).find_first_not_of(" \t\r\n") != std::string::npos)) {
      error = "Text is not permitted in this element.";
      return false;
    }
  }
  return true;
}

xmlNodePtr Constraints::create(const std::string &name,
                               std::set<xmlAttrPtr> &pending) const {

  // Find a finite, minimal derivation by relaxation. Recursive declarations
  // with no terminating alternative never acquire a derivation.
  std::map<std::string, std::vector<std::string>> plans;
  std::map<std::string, std::size_t> costs;

  constexpr std::size_t limit = 4096;

  struct Plan {
    bool ok = false;
    std::size_t cost = 0;
    std::vector<std::string> children;
  };

  std::function<Plan(xmlElementContentPtr)> derive = [&](auto c) -> Plan {
    if (!c || c->ocur == XML_ELEMENT_CONTENT_OPT ||
        c->ocur == XML_ELEMENT_CONTENT_MULT)
      return {true, 0, {}};
    if (c->type == XML_ELEMENT_CONTENT_PCDATA)
      return {true, 0, {}};
    if (c->type == XML_ELEMENT_CONTENT_ELEMENT) {
      std::string key(reinterpret_cast<const char *>(c->name));
      auto i = costs.find(key);
      return i == costs.end() ? Plan{} : Plan{true, i->second, {key}};
    }
    auto left = derive(c->c1), right = derive(c->c2);
    if (c->type == XML_ELEMENT_CONTENT_OR)
      return !left.ok ? right
                      : (!right.ok || left.cost <= right.cost ? left : right);
    if (!left.ok || !right.ok || left.cost + right.cost >= limit)
      return {};
    left.cost += right.cost;
    left.children.insert(left.children.end(), right.children.begin(),
                         right.children.end());
    return left;
  };

  auto declared = names();
  bool changed = true;

  while (changed) {
    changed = false;
    for (const auto &key : declared) {
      auto e = element(BAD_CAST key.c_str());
      auto p = (e->etype == XML_ELEMENT_TYPE_EMPTY ||
                e->etype == XML_ELEMENT_TYPE_ANY)
                   ? Plan{true, 0, {}}
                   : derive(e->content);
      if (p.ok && p.cost < limit &&
          (!costs.count(key) || p.cost + 1 < costs[key])) {
        costs[key] = p.cost + 1;
        plans[key] = p.children;
        changed = true;
      }
    }
  }

  if (!plans.count(name))
    return nullptr;

  std::function<xmlNodePtr(const std::string &)> build = [&](const auto &key) {
    auto n = xmlNewDocNode(doc, nullptr, BAD_CAST key.c_str(), nullptr);
    for (auto a : attributes(n->name)) {
      if (a->defaultValue)
        xmlSetProp(n, a->name, a->defaultValue);
      else if (a->def == XML_ATTRIBUTE_REQUIRED)
        pending.insert(xmlSetProp(n, a->name, BAD_CAST ""));
    }
    for (const auto &child : plans[key])
      xmlAddChild(n, build(child));
    return n;
  };

  return build(name);
}
} // namespace editor
