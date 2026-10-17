#include "elff/io/C/xmlDTDValidator.h"

#include "elff/config/_config.hpp"

#include <limits.h>
#include <stdio.h>
#include <string.h>

#ifdef ELFF_USE_LIBXML2
#include <libxml/parser.h>
#include <libxml/tree.h>
#include <libxml/xmlerror.h>
#endif

typedef struct
{
  char* buffer;
  size_t capacity;
} elff_xml_diagnostic;

static void
elff_xml_set_diagnostic(elff_xml_diagnostic* diagnostic, const char* message)
{
  size_t message_length;

  if (diagnostic->buffer == NULL || diagnostic->capacity == 0) {
    return;
  }

  if (message == NULL) {
    diagnostic->buffer[0] = '\0';
    return;
  }

  message_length = strlen(message);
  while (message_length > 0 &&
         (message[message_length - 1] == '\n' ||
          message[message_length - 1] == '\r')) {
    --message_length;
  }

  (void)snprintf(diagnostic->buffer,
                 diagnostic->capacity,
                 "%.*s",
                 (int)(message_length > INT_MAX ? INT_MAX : message_length),
                 message);
}

#ifdef ELFF_USE_LIBXML2
static void
elff_xml_capture_error(void* user_data, const xmlError* error)
{
  elff_xml_diagnostic* diagnostic = (elff_xml_diagnostic*)user_data;

  if (error == NULL || error->level < XML_ERR_ERROR ||
      diagnostic->buffer == NULL || diagnostic->capacity == 0 ||
      diagnostic->buffer[0] != '\0') {
    return;
  }

  elff_xml_set_diagnostic(diagnostic, error->message);
}

static int
elff_xml_diagnostic_is_empty(const elff_xml_diagnostic* diagnostic)
{
  return diagnostic->buffer == NULL || diagnostic->capacity == 0 ||
         diagnostic->buffer[0] == '\0';
}

static int
elff_xml_node_has_unqualified_name(const xmlNode* node, const char* name)
{
  return node != NULL && node->type == XML_ELEMENT_NODE && node->ns == NULL &&
         xmlStrEqual(node->name, (const xmlChar*)name);
}

static elff_xml_dtd_status
elff_xml_validate_document(xmlDoc* document,
                           const char* dtd,
                           size_t dtd_size,
                           elff_xml_diagnostic* diagnostic)
{
  xmlParserCtxt* dtd_context = NULL;
  xmlParserCtxt* validation_context = NULL;
  xmlParserInput* dtd_input = NULL;
  xmlDtd* parsed_dtd = NULL;
  elff_xml_dtd_status status = ELFF_XML_DTD_INTERNAL_ERROR;
  int parser_options = XML_PARSE_NONET | XML_PARSE_NO_XXE;

  dtd_context = xmlNewParserCtxt();
  validation_context = xmlNewParserCtxt();
  if (dtd_context == NULL || validation_context == NULL) {
    elff_xml_set_diagnostic(diagnostic,
                            "Unable to allocate a libxml2 parser context.");
    goto cleanup;
  }

  xmlCtxtSetErrorHandler(dtd_context, elff_xml_capture_error, diagnostic);
  (void)xmlCtxtSetOptions(dtd_context, parser_options | XML_PARSE_DTDLOAD);
  dtd_input = xmlNewInputFromMemory(
    "in-memory.dtd", dtd, dtd_size, XML_INPUT_BUF_STATIC);
  if (dtd_input == NULL) {
    elff_xml_set_diagnostic(diagnostic,
                            "Unable to allocate a libxml2 DTD input.");
    goto cleanup;
  }

  parsed_dtd = xmlCtxtParseDtd(dtd_context,
                               dtd_input,
                               NULL,
                               (const xmlChar*)"in-memory.dtd");
  if (parsed_dtd == NULL) {
    if (elff_xml_diagnostic_is_empty(diagnostic)) {
      elff_xml_set_diagnostic(diagnostic, "Unable to parse the DTD.");
    }
    status = ELFF_XML_DTD_MALFORMED_DTD;
    goto cleanup;
  }

  xmlCtxtSetErrorHandler(
    validation_context, elff_xml_capture_error, diagnostic);
  if (xmlCtxtValidateDtd(validation_context, document, parsed_dtd) == 1) {
    elff_xml_set_diagnostic(diagnostic, NULL);
    status = ELFF_XML_DTD_VALID;
  } else {
    if (elff_xml_diagnostic_is_empty(diagnostic)) {
      elff_xml_set_diagnostic(
        diagnostic, "The XML document does not conform to the DTD.");
    }
    status = ELFF_XML_DTD_INVALID_DOCUMENT;
  }

cleanup:
  xmlFreeDtd(parsed_dtd);
  xmlFreeParserCtxt(validation_context);
  xmlFreeParserCtxt(dtd_context);
  return status;
}
#endif

elff_xml_dtd_status
elff_xml_dtd_validate_memory(const char* xml,
                             size_t xml_size,
                             const char* dtd,
                             size_t dtd_size,
                             char* diagnostic_buffer,
                             size_t diagnostic_capacity)
{
  elff_xml_diagnostic diagnostic = { diagnostic_buffer,
                                     diagnostic_capacity };

  elff_xml_set_diagnostic(&diagnostic, NULL);

  if (xml == NULL || dtd == NULL || xml_size > INT_MAX) {
    elff_xml_set_diagnostic(&diagnostic, "Invalid XML DTD validator argument.");
    return ELFF_XML_DTD_INVALID_ARGUMENT;
  }

#ifndef ELFF_USE_LIBXML2
  (void)xml_size;
  (void)dtd_size;
  elff_xml_set_diagnostic(&diagnostic,
                          "ELFF was built without libxml2 support.");
  return ELFF_XML_DTD_UNAVAILABLE;
#else
  xmlParserCtxt* xml_context = NULL;
  xmlDoc* document = NULL;
  elff_xml_dtd_status status = ELFF_XML_DTD_INTERNAL_ERROR;
  int parser_options = XML_PARSE_NONET | XML_PARSE_NO_XXE;

  xml_context = xmlNewParserCtxt();
  if (xml_context == NULL) {
    elff_xml_set_diagnostic(&diagnostic,
                            "Unable to allocate a libxml2 parser context.");
    goto cleanup;
  }

  xmlCtxtSetErrorHandler(xml_context, elff_xml_capture_error, &diagnostic);
  document = xmlCtxtReadMemory(xml_context,
                               xml,
                               (int)xml_size,
                               "in-memory.xml",
                               NULL,
                               parser_options);
  if (document == NULL) {
    if (elff_xml_diagnostic_is_empty(&diagnostic)) {
      elff_xml_set_diagnostic(&diagnostic, "Unable to parse the XML document.");
    }
    status = ELFF_XML_DTD_MALFORMED_XML;
    goto cleanup;
  }

  status =
    elff_xml_validate_document(document, dtd, dtd_size, &diagnostic);

cleanup:
  xmlFreeDoc(document);
  xmlFreeParserCtxt(xml_context);
  return status;
#endif
}

elff_xml_dtd_status
elff_xml_dtd_validate_node(elff_xml_node_t node,
                           const char* dtd,
                           size_t dtd_size,
                           char* diagnostic_buffer,
                           size_t diagnostic_capacity)
{
  elff_xml_diagnostic diagnostic = { diagnostic_buffer,
                                     diagnostic_capacity };

  elff_xml_set_diagnostic(&diagnostic, NULL);

  if (node == NULL || dtd == NULL) {
    elff_xml_set_diagnostic(&diagnostic, "Invalid XML DTD validator argument.");
    return ELFF_XML_DTD_INVALID_ARGUMENT;
  }

#ifndef ELFF_USE_LIBXML2
  (void)dtd_size;
  elff_xml_set_diagnostic(&diagnostic,
                          "ELFF was built without libxml2 support.");
  return ELFF_XML_DTD_UNAVAILABLE;
#else
  const xmlNode* selected_element = (const xmlNode*)node;
  xmlDoc* subtree_document = NULL;
  xmlNode* selected_copy = NULL;
  elff_xml_dtd_status status = ELFF_XML_DTD_INTERNAL_ERROR;

  if (selected_element->type != XML_ELEMENT_NODE) {
    elff_xml_set_diagnostic(&diagnostic,
                            "The selected XML node is not an element.");
    return ELFF_XML_DTD_INVALID_ARGUMENT;
  }

  subtree_document = xmlNewDoc((const xmlChar*)"1.0");
  if (subtree_document == NULL) {
    elff_xml_set_diagnostic(
      &diagnostic, "Unable to allocate a document for the selected subtree.");
    goto cleanup;
  }

  selected_copy = xmlDocCopyNode((xmlNode*)selected_element,
                                 subtree_document,
                                 1);
  if (selected_copy == NULL) {
    elff_xml_set_diagnostic(
      &diagnostic, "Unable to copy the selected XML subtree.");
    goto cleanup;
  }

  (void)xmlDocSetRootElement(subtree_document, selected_copy);
  if (xmlReconciliateNs(subtree_document, selected_copy) < 0) {
    elff_xml_set_diagnostic(
      &diagnostic, "Unable to reconcile namespaces in the selected subtree.");
    goto cleanup;
  }

  status =
    elff_xml_validate_document(subtree_document, dtd, dtd_size, &diagnostic);

cleanup:
  xmlFreeDoc(subtree_document);
  return status;
#endif
}

elff_xml_dtd_status
elff_xml_dtd_validate_subtree_memory(const char* xml,
                                     size_t xml_size,
                                     const char* container_element_name,
                                     const char* subtree_element_name,
                                     const char* dtd,
                                     size_t dtd_size,
                                     char* diagnostic_buffer,
                                     size_t diagnostic_capacity)
{
  elff_xml_diagnostic diagnostic = { diagnostic_buffer,
                                     diagnostic_capacity };

  elff_xml_set_diagnostic(&diagnostic, NULL);

  if (xml == NULL || container_element_name == NULL ||
      subtree_element_name == NULL || dtd == NULL || xml_size > INT_MAX ||
      container_element_name[0] == '\0' || subtree_element_name[0] == '\0') {
    elff_xml_set_diagnostic(&diagnostic, "Invalid XML DTD validator argument.");
    return ELFF_XML_DTD_INVALID_ARGUMENT;
  }

#ifndef ELFF_USE_LIBXML2
  (void)xml_size;
  (void)dtd_size;
  elff_xml_set_diagnostic(&diagnostic,
                          "ELFF was built without libxml2 support.");
  return ELFF_XML_DTD_UNAVAILABLE;
#else
  xmlParserCtxt* xml_context = NULL;
  xmlDoc* document = NULL;
  xmlDoc* subtree_document = NULL;
  xmlNode* document_element = NULL;
  xmlNode* selected_element = NULL;
  xmlNode* selected_copy = NULL;
  xmlNode* child = NULL;
  size_t matching_children = 0;
  elff_xml_dtd_status status = ELFF_XML_DTD_INTERNAL_ERROR;
  int parser_options = XML_PARSE_NONET | XML_PARSE_NO_XXE;

  xml_context = xmlNewParserCtxt();
  if (xml_context == NULL) {
    elff_xml_set_diagnostic(&diagnostic,
                            "Unable to allocate a libxml2 parser context.");
    goto cleanup;
  }

  xmlCtxtSetErrorHandler(xml_context, elff_xml_capture_error, &diagnostic);
  document = xmlCtxtReadMemory(xml_context,
                               xml,
                               (int)xml_size,
                               "in-memory.xml",
                               NULL,
                               parser_options);
  if (document == NULL) {
    if (elff_xml_diagnostic_is_empty(&diagnostic)) {
      elff_xml_set_diagnostic(&diagnostic, "Unable to parse the XML document.");
    }
    status = ELFF_XML_DTD_MALFORMED_XML;
    goto cleanup;
  }

  document_element = xmlDocGetRootElement(document);
  if (!elff_xml_node_has_unqualified_name(document_element,
                                           container_element_name)) {
    if (diagnostic.buffer != NULL && diagnostic.capacity > 0) {
      (void)snprintf(diagnostic.buffer,
                     diagnostic.capacity,
                     "Expected <%s> as the document element.",
                     container_element_name);
    }
    status = ELFF_XML_DTD_INVALID_DOCUMENT;
    goto cleanup;
  }

  for (child = document_element->children; child != NULL; child = child->next) {
    if (elff_xml_node_has_unqualified_name(child, subtree_element_name)) {
      selected_element = child;
      ++matching_children;
    }
  }

  if (matching_children == 0) {
    if (diagnostic.buffer != NULL && diagnostic.capacity > 0) {
      (void)snprintf(diagnostic.buffer,
                     diagnostic.capacity,
                     "The <%s> element has no direct <%s> child.",
                     container_element_name,
                     subtree_element_name);
    }
    status = ELFF_XML_DTD_SUBTREE_NOT_FOUND;
    goto cleanup;
  }

  if (matching_children != 1) {
    if (diagnostic.buffer != NULL && diagnostic.capacity > 0) {
      (void)snprintf(diagnostic.buffer,
                     diagnostic.capacity,
                     "The <%s> element has more than one direct <%s> child.",
                     container_element_name,
                     subtree_element_name);
    }
    status = ELFF_XML_DTD_SUBTREE_NOT_UNIQUE;
    goto cleanup;
  }

  subtree_document = xmlNewDoc((const xmlChar*)"1.0");
  if (subtree_document == NULL) {
    elff_xml_set_diagnostic(
      &diagnostic, "Unable to allocate a document for the selected subtree.");
    goto cleanup;
  }

  selected_copy = xmlDocCopyNode(selected_element, subtree_document, 1);
  if (selected_copy == NULL) {
    elff_xml_set_diagnostic(
      &diagnostic, "Unable to copy the selected XML subtree.");
    goto cleanup;
  }

  (void)xmlDocSetRootElement(subtree_document, selected_copy);
  if (xmlReconciliateNs(subtree_document, selected_copy) < 0) {
    elff_xml_set_diagnostic(
      &diagnostic, "Unable to reconcile namespaces in the selected subtree.");
    goto cleanup;
  }

  status =
    elff_xml_validate_document(subtree_document, dtd, dtd_size, &diagnostic);

cleanup:
  xmlFreeDoc(subtree_document);
  xmlFreeDoc(document);
  xmlFreeParserCtxt(xml_context);
  return status;
#endif
}

const char*
elff_xml_dtd_status_string(elff_xml_dtd_status status)
{
  switch (status) {
    case ELFF_XML_DTD_VALID:
      return "valid";
    case ELFF_XML_DTD_INVALID_DOCUMENT:
      return "invalid document";
    case ELFF_XML_DTD_MALFORMED_XML:
      return "malformed XML";
    case ELFF_XML_DTD_MALFORMED_DTD:
      return "malformed DTD";
    case ELFF_XML_DTD_INVALID_ARGUMENT:
      return "invalid argument";
    case ELFF_XML_DTD_INTERNAL_ERROR:
      return "internal error";
    case ELFF_XML_DTD_UNAVAILABLE:
      return "libxml2 unavailable";
    case ELFF_XML_DTD_SUBTREE_NOT_FOUND:
      return "subtree not found";
    case ELFF_XML_DTD_SUBTREE_NOT_UNIQUE:
      return "subtree not unique";
    default:
      return "unknown status";
  }
}
