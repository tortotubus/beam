#include "elff/config/_config.hpp"
#include "elff/io/C/xmlCaseParameters.h"
#ifdef ELFF_USE_LIBXML2
#include "../C/xmlCaseInternal.h"
#include <libxml/xmlsave.h>
#endif
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <limits>
#include <locale>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace {
using Document = std::unique_ptr<elff_xml_document_handle,
                                 decltype(&elff_xml_document_destroy)>;
using Parameters = std::unique_ptr<elff_xml_parameters_handle,
                                   decltype(&elff_xml_parameters_destroy)>;
#ifdef ELFF_USE_LIBXML2
elff_xml_status render_values(xmlNode *root,
                                  elff_xml_parameters_t parameters,
                                  elff_xml_error *error) {
  std::vector<xmlNode *> pending{root};
  while (!pending.empty()) {
    xmlNode *node = pending.back();
    pending.pop_back();
    // These sections are source definitions, not downstream attributes.
    if (node->parent == root &&
        (xmlStrEqual(node->name, BAD_CAST "Parameters") ||
         xmlStrEqual(node->name, BAD_CAST "Nondimensionalization")))
      continue;
    for (xmlNode *child = node->children; child; child = child->next)
      if (child->type == XML_ELEMENT_NODE)
        pending.push_back(child);
    std::vector<xmlAttr*> slots;
    bool leaf = true;
    for (xmlNode* child = node->children; child; child = child->next)
      if (child->type == XML_ELEMENT_NODE) leaf = false;
    if (leaf) slots.push_back(nullptr);
    for (xmlAttr* attr = node->properties; attr; attr = attr->next)
      slots.push_back(attr);
    for (xmlAttr* attr : slots) {
      std::unique_ptr<xmlChar, decltype(xmlFree)> content(
          attr ? xmlNodeListGetString(node->doc, attr->children, 1) : xmlNodeGetContent(node), xmlFree);
      if (!content)
        continue;
      std::string source(reinterpret_cast<const char *>(content.get()));
      if (source.find("${") == std::string::npos)
        continue;
      const auto first = source.find_first_not_of(" \t\r\n");
      const auto last = source.find_last_not_of(" \t\r\n");
      std::string expression;
      elff_xml_error detail{};
      elff_xml_status status = ELFF_XML_INVALID_VALUE;
      double value = 0;
      if (source.size() > 4099 || source.compare(first, 2, "${") != 0 ||
          source[last] != '}' || source.find('}', first) != last) {
        elff_xml_error_set(&detail, 0, 0,
                           "Expected a whole-value ${expression}, without "
                           "string interpolation.");
      } else if (!parameters) {
        elff_xml_error_set(&detail, 0, 0,
                           "Render expressions require a Parameters section.");
      } else {
        expression = source.substr(first + 2, last - first - 2);
        status = elff_xml_parameters_evaluate(parameters, expression.c_str(),
                                              &value, &detail);
      }
      if (status != ELFF_XML_OK) {
        std::string message =
            "<" + std::string(reinterpret_cast<const char *>(node->name)) +
            (attr ? "> attribute '" + std::string(reinterpret_cast<const char *>(attr->name)) + "'"
                  : "> text") +
            " (" + source.substr(0, 128) + "): " + detail.message;
        elff_xml_error_set(error, static_cast<size_t>(xmlGetLineNo(node)), 0,
                           message.c_str());
        return status;
      }
      std::ostringstream number;
      number.imbue(std::locale::classic());
      number << std::setprecision(std::numeric_limits<double>::max_digits10)
             << value;
      if (!attr) {
        xmlNodeSetContent(node, BAD_CAST number.str().c_str());
      } else if (!xmlSetNsProp(node, attr->ns, attr->name,
                        BAD_CAST number.str().c_str())) {
        elff_xml_error_set(error, xmlGetLineNo(node), 0,
                           "Unable to replace rendered attribute.");
        return ELFF_XML_INTERNAL_ERROR;
      }
    }
  }
  return ELFF_XML_OK;
}
#endif
} // namespace

extern "C" elff_xml_status
elff_xml_case_render(elff_xml_document_t source,
                     elff_xml_document_t *out_document, elff_xml_error *error) {
  elff_xml_error_clear(error);
  if (out_document)
    *out_document = nullptr;
  if (!source || !out_document) {
    elff_xml_error_set(error, 0, 0, "Source document and output are required.");
    return ELFF_XML_INVALID_ARGUMENT;
  }
#ifndef ELFF_USE_LIBXML2
  elff_xml_error_set(error, 0, 0, "XML support is unavailable in this build.");
  return ELFF_XML_UNAVAILABLE;
#else
  try {
    elff_xml_node_t root = nullptr, section = nullptr;
    auto status = elff_xml_document_root(source, "Case", &root, error);
    if (status != ELFF_XML_OK)
      return status;
    int version = 0;
    elff_xml_error version_error{};
    const bool v2 = elff_xml_node_attribute_int(root, "formatVersion", &version,
                                               &version_error) == ELFF_XML_OK && version == 2;
    status = elff_xml_node_unique_child(root, "Parameters", &section, error);
    Parameters parameters(nullptr, elff_xml_parameters_destroy);
    if (status == ELFF_XML_OK || (v2 && status == ELFF_XML_NOT_FOUND)) {
      elff_xml_parameters_t context = nullptr;
      status = elff_xml_parameters_create(root, &context, error);
      if (status != ELFF_XML_OK)
        return status;
      parameters.reset(context);
    } else if (status == ELFF_XML_NOT_FOUND) {
      status = elff_xml_node_unique_child(root, "Nondimensionalization",
                                          &section, error);
      if (status == ELFF_XML_OK) {
        elff_xml_error_set(error, 0, 0,
                           "Nondimensionalization requires Parameters.");
        return ELFF_XML_INVALID_VALUE;
      }
      if (status != ELFF_XML_NOT_FOUND)
        return status;
      elff_xml_error_clear(error);
    } else
      return status;
    Document rendered(
        static_cast<elff_xml_document_t>(std::calloc(1, sizeof(*source))),
        elff_xml_document_destroy);
    if (!rendered)
      throw std::bad_alloc();
    rendered->document = xmlCopyDoc(source->document, 1);
    if (!rendered->document)
      throw std::bad_alloc();
    status = render_values(xmlDocGetRootElement(rendered->document),
                               parameters.get(), error);
    if (status != ELFF_XML_OK)
      return status;
    *out_document = rendered.release();
    return ELFF_XML_OK;
  } catch (...) {
    elff_xml_error_set(error, 0, 0, "Unable to allocate or render case XML.");
    return ELFF_XML_INTERNAL_ERROR;
  }
#endif
}

extern "C" elff_xml_status
elff_xml_document_write_file(elff_xml_document_t document, const char *path,
                             elff_xml_error *error) {
  elff_xml_error_clear(error);
  if (!document || !path || !*path) {
    elff_xml_error_set(error, 0, 0,
                       "Document and nonempty output path are required.");
    return ELFF_XML_INVALID_ARGUMENT;
  }
#ifndef ELFF_USE_LIBXML2
  elff_xml_error_set(error, 0, 0, "XML support is unavailable in this build.");
  return ELFF_XML_UNAVAILABLE;
#else
  if (xmlSaveFormatFileEnc(path, document->document, "UTF-8", 1) < 0) {
    elff_xml_error_set(error, 0, 0, "Unable to write rendered XML diagnostic.");
    return ELFF_XML_IO_ERROR;
  }
  return ELFF_XML_OK;
#endif
}

extern "C" elff_xml_status
elff_xml_case_load_rendered_file(const char *path, int write_diagnostic,
                                 elff_xml_document_t *out_document,
                                 elff_xml_error *error) {
  elff_xml_error_clear(error);
  if (out_document)
    *out_document = nullptr;
  if (!path || !*path || !out_document) {
    elff_xml_error_set(error, 0, 0,
                       "Source path and output document are required.");
    return ELFF_XML_INVALID_ARGUMENT;
  }
  try {
    elff_xml_document_t source = nullptr;
    auto status = elff_xml_document_load_file(path, &source, error);
    if (status != ELFF_XML_OK)
      return status;
    Document owner(source, elff_xml_document_destroy);
    status = elff_xml_case_render(source, out_document, error);
    if (status != ELFF_XML_OK)
      return status;
  } catch (...) {
    elff_xml_error_set(error, 0, 0, "Unable to load and render case XML.");
    return ELFF_XML_INTERNAL_ERROR;
  }
  if (write_diagnostic) {
    // Diagnostics are strictly best effort, including path allocation failures.
    try {
      std::string output(path);
      if (output.size() >= 4 &&
          output.compare(output.size() - 4, 4, ".xml") == 0)
        output.resize(output.size() - 4);
      output += ".rendered.xml";
      elff_xml_error diagnostic{};
      if (elff_xml_document_write_file(*out_document, output.c_str(),
                                       &diagnostic) != ELFF_XML_OK)
        std::fprintf(stderr, "Warning: %s (%s)\n", diagnostic.message,
                     output.c_str());
    } catch (...) {
      std::fprintf(stderr,
                   "Warning: unable to write rendered XML diagnostic for %s.\n",
                   path);
    }
  }
  return ELFF_XML_OK;
}
