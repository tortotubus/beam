#include "elff/io/C/xmlCase.h"

#include "elff/config/_config.hpp"

#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef ELFF_USE_LIBXML2
#include <libxml/parser.h>
#include <libxml/tree.h>
#include <libxml/xmlerror.h>

struct elff_xml_document_handle
{
  xmlDoc* document;
};

typedef struct
{
  elff_xml_error* error;
  int saw_io_error;
} elff_xml_error_capture;
#endif

void
elff_xml_error_clear(elff_xml_error* error)
{
  if (error == NULL) {
    return;
  }

  error->line = 0;
  error->column = 0;
  error->message[0] = '\0';
}

void
elff_xml_error_set(elff_xml_error* error,
                   size_t line,
                   size_t column,
                   const char* message)
{
  size_t length;

  if (error == NULL) {
    return;
  }

  error->line = line;
  error->column = column;
  if (message == NULL) {
    error->message[0] = '\0';
    return;
  }

  length = strlen(message);
  while (length > 0 &&
         (message[length - 1] == '\n' || message[length - 1] == '\r')) {
    --length;
  }

  (void)snprintf(error->message,
                 sizeof(error->message),
                 "%.*s",
                 (int)(length > INT_MAX ? INT_MAX : length),
                 message);
}

#ifdef ELFF_USE_LIBXML2
static void
elff_xml_capture_error(void* user_data, const xmlError* xml_error)
{
  elff_xml_error_capture* capture = (elff_xml_error_capture*)user_data;

  if (capture == NULL || xml_error == NULL ||
      xml_error->level < XML_ERR_ERROR) {
    return;
  }

  if (xml_error->domain == XML_FROM_IO) {
    capture->saw_io_error = 1;
  }

  if (capture->error == NULL || capture->error->message[0] != '\0') {
    return;
  }

  elff_xml_error_set(capture->error,
                     xml_error->line > 0 ? (size_t)xml_error->line : 0,
                     xml_error->int2 > 0 ? (size_t)xml_error->int2 : 0,
                     xml_error->message);
}

static const xmlNode*
elff_xml_node_from_handle(elff_xml_node_t node)
{
  return (const xmlNode*)node;
}

static elff_xml_node_t
elff_xml_node_to_handle(const xmlNode* node)
{
  return (elff_xml_node_t)node;
}

static size_t
elff_xml_node_line_number(const xmlNode* node)
{
  long line = xmlGetLineNo(node);
  return line > 0 ? (size_t)line : 0;
}

static int
elff_xml_node_has_unqualified_name(const xmlNode* node, const char* name)
{
  return node != NULL && node->type == XML_ELEMENT_NODE && node->ns == NULL &&
         xmlStrEqual(node->name, (const xmlChar*)name);
}

static elff_xml_status
elff_xml_document_wrap(xmlDoc* parsed,
                       elff_xml_document_t* out_document,
                       elff_xml_error* error)
{
  elff_xml_document_t document;

  document = (elff_xml_document_t)malloc(sizeof(*document));
  if (document == NULL) {
    xmlFreeDoc(parsed);
    elff_xml_error_set(error, 0, 0, "Unable to allocate an XML document handle.");
    return ELFF_XML_INTERNAL_ERROR;
  }

  document->document = parsed;
  *out_document = document;
  return ELFF_XML_OK;
}

static const char*
elff_xml_skip_space(const char* value)
{
  while (*value != '\0' && isspace((unsigned char)*value)) {
    ++value;
  }
  return value;
}

static int
elff_xml_only_trailing_space(const char* value)
{
  return *elff_xml_skip_space(value) == '\0';
}
#endif

elff_xml_status
elff_xml_document_load_memory(const char* xml,
                              size_t xml_size,
                              elff_xml_document_t* out_document,
                              elff_xml_error* error)
{
  elff_xml_error_clear(error);

  if (out_document != NULL) {
    *out_document = NULL;
  }
  if (xml == NULL || out_document == NULL || xml_size > INT_MAX) {
    elff_xml_error_set(error, 0, 0, "Invalid XML document argument.");
    return ELFF_XML_INVALID_ARGUMENT;
  }

#ifndef ELFF_USE_LIBXML2
  (void)xml_size;
  elff_xml_error_set(error, 0, 0, "ELFF was built without libxml2 support.");
  return ELFF_XML_UNAVAILABLE;
#else
  xmlParserCtxt* context = NULL;
  xmlDoc* parsed = NULL;
  elff_xml_error_capture capture = { error, 0 };
  int options = XML_PARSE_NONET | XML_PARSE_NO_XXE;

  context = xmlNewParserCtxt();
  if (context == NULL) {
    elff_xml_error_set(error, 0, 0, "Unable to allocate a libxml2 parser context.");
    return ELFF_XML_INTERNAL_ERROR;
  }

  xmlCtxtSetErrorHandler(context, elff_xml_capture_error, &capture);
  parsed = xmlCtxtReadMemory(
    context, xml, (int)xml_size, "in-memory.xml", NULL, options);
  xmlFreeParserCtxt(context);

  if (parsed == NULL) {
    if (error == NULL || error->message[0] == '\0') {
      elff_xml_error_set(error, 0, 0, "Unable to parse the XML document.");
    }
    return ELFF_XML_MALFORMED_DOCUMENT;
  }

  return elff_xml_document_wrap(parsed, out_document, error);
#endif
}

elff_xml_status
elff_xml_document_load_file(const char* path,
                            elff_xml_document_t* out_document,
                            elff_xml_error* error)
{
  elff_xml_error_clear(error);

  if (out_document != NULL) {
    *out_document = NULL;
  }
  if (path == NULL || path[0] == '\0' || out_document == NULL) {
    elff_xml_error_set(error, 0, 0, "Invalid XML file argument.");
    return ELFF_XML_INVALID_ARGUMENT;
  }

#ifndef ELFF_USE_LIBXML2
  elff_xml_error_set(error, 0, 0, "ELFF was built without libxml2 support.");
  return ELFF_XML_UNAVAILABLE;
#else
  FILE* input = NULL;
  xmlParserCtxt* context = NULL;
  xmlDoc* parsed = NULL;
  elff_xml_error_capture capture = { error, 0 };
  int options = XML_PARSE_NONET | XML_PARSE_NO_XXE;

  input = fopen(path, "rb");
  if (input == NULL) {
    int open_error = errno;
    char message[512];

    (void)snprintf(message,
                   sizeof(message),
                   "Unable to open XML file '%s': %s.",
                   path,
                   strerror(open_error));
    elff_xml_error_set(error, 0, 0, message);
    return ELFF_XML_IO_ERROR;
  }
  (void)fclose(input);

  context = xmlNewParserCtxt();
  if (context == NULL) {
    elff_xml_error_set(error, 0, 0, "Unable to allocate a libxml2 parser context.");
    return ELFF_XML_INTERNAL_ERROR;
  }

  xmlCtxtSetErrorHandler(context, elff_xml_capture_error, &capture);
  parsed = xmlCtxtReadFile(context, path, NULL, options);
  xmlFreeParserCtxt(context);

  if (parsed == NULL) {
    if (error == NULL || error->message[0] == '\0') {
      elff_xml_error_set(error, 0, 0, "Unable to parse the XML file.");
    }
    return capture.saw_io_error ? ELFF_XML_IO_ERROR
                                : ELFF_XML_MALFORMED_DOCUMENT;
  }

  return elff_xml_document_wrap(parsed, out_document, error);
#endif
}

void
elff_xml_document_destroy(elff_xml_document_t document)
{
#ifdef ELFF_USE_LIBXML2
  if (document != NULL) {
    xmlFreeDoc(document->document);
    free(document);
  }
#else
  (void)document;
#endif
}

elff_xml_status
elff_xml_document_root(elff_xml_document_t document,
                       const char* expected_name,
                       elff_xml_node_t* out_root,
                       elff_xml_error* error)
{
  elff_xml_error_clear(error);

  if (out_root != NULL) {
    *out_root = NULL;
  }
  if (document == NULL || expected_name == NULL || expected_name[0] == '\0' ||
      out_root == NULL) {
    elff_xml_error_set(error, 0, 0, "Invalid XML root-element argument.");
    return ELFF_XML_INVALID_ARGUMENT;
  }

#ifndef ELFF_USE_LIBXML2
  elff_xml_error_set(error, 0, 0, "ELFF was built without libxml2 support.");
  return ELFF_XML_UNAVAILABLE;
#else
  const xmlNode* root = xmlDocGetRootElement(document->document);

  if (root == NULL) {
    elff_xml_error_set(error, 0, 0, "The XML document has no document element.");
    return ELFF_XML_NOT_FOUND;
  }
  if (!elff_xml_node_has_unqualified_name(root, expected_name)) {
    char message[ELFF_XML_ERROR_MESSAGE_CAPACITY];
    (void)snprintf(message,
                   sizeof(message),
                   "Expected <%s> as the document element.",
                   expected_name);
    elff_xml_error_set(error, elff_xml_node_line_number(root), 0, message);
    return ELFF_XML_INVALID_VALUE;
  }

  *out_root = elff_xml_node_to_handle(root);
  return ELFF_XML_OK;
#endif
}

elff_xml_status
elff_xml_node_unique_child(elff_xml_node_t parent,
                           const char* child_name,
                           elff_xml_node_t* out_child,
                           elff_xml_error* error)
{
  elff_xml_error_clear(error);

  if (out_child != NULL) {
    *out_child = NULL;
  }
  if (parent == NULL || child_name == NULL || child_name[0] == '\0' ||
      out_child == NULL) {
    elff_xml_error_set(error, 0, 0, "Invalid XML child-selection argument.");
    return ELFF_XML_INVALID_ARGUMENT;
  }

#ifndef ELFF_USE_LIBXML2
  elff_xml_error_set(error, 0, 0, "ELFF was built without libxml2 support.");
  return ELFF_XML_UNAVAILABLE;
#else
  const xmlNode* parent_node = elff_xml_node_from_handle(parent);
  const xmlNode* child;
  const xmlNode* selected = NULL;
  size_t count = 0;

  for (child = parent_node->children; child != NULL; child = child->next) {
    if (elff_xml_node_has_unqualified_name(child, child_name)) {
      selected = child;
      ++count;
    }
  }

  if (count == 0) {
    char message[ELFF_XML_ERROR_MESSAGE_CAPACITY];
    (void)snprintf(message,
                   sizeof(message),
                   "The selected element has no direct <%s> child.",
                   child_name);
    elff_xml_error_set(
      error, elff_xml_node_line_number(parent_node), 0, message);
    return ELFF_XML_NOT_FOUND;
  }
  if (count != 1) {
    char message[ELFF_XML_ERROR_MESSAGE_CAPACITY];
    (void)snprintf(message,
                   sizeof(message),
                   "The selected element has more than one direct <%s> child.",
                   child_name);
    elff_xml_error_set(
      error, elff_xml_node_line_number(parent_node), 0, message);
    return ELFF_XML_NOT_UNIQUE;
  }

  *out_child = elff_xml_node_to_handle(selected);
  return ELFF_XML_OK;
#endif
}

elff_xml_status
elff_xml_node_element_child_count(elff_xml_node_t parent,
                                  size_t* out_count,
                                  elff_xml_error* error)
{
  elff_xml_error_clear(error);
  if (out_count != NULL) {
    *out_count = 0;
  }
  if (parent == NULL || out_count == NULL) {
    elff_xml_error_set(error, 0, 0, "Invalid XML child-count argument.");
    return ELFF_XML_INVALID_ARGUMENT;
  }
#ifndef ELFF_USE_LIBXML2
  elff_xml_error_set(error, 0, 0, "ELFF was built without libxml2 support.");
  return ELFF_XML_UNAVAILABLE;
#else
  const xmlNode* child;
  for (child = elff_xml_node_from_handle(parent)->children; child != NULL;
       child = child->next) {
    if (child->type == XML_ELEMENT_NODE) {
      ++*out_count;
    }
  }
  return ELFF_XML_OK;
#endif
}

elff_xml_status
elff_xml_node_element_child_at(elff_xml_node_t parent,
                               size_t index,
                               elff_xml_node_t* out_child,
                               elff_xml_error* error)
{
  elff_xml_error_clear(error);
  if (out_child != NULL) {
    *out_child = NULL;
  }
  if (parent == NULL || out_child == NULL) {
    elff_xml_error_set(error, 0, 0, "Invalid XML child-index argument.");
    return ELFF_XML_INVALID_ARGUMENT;
  }
#ifndef ELFF_USE_LIBXML2
  (void)index;
  elff_xml_error_set(error, 0, 0, "ELFF was built without libxml2 support.");
  return ELFF_XML_UNAVAILABLE;
#else
  const xmlNode* parent_node = elff_xml_node_from_handle(parent);
  const xmlNode* child;
  size_t current = 0;
  for (child = parent_node->children; child != NULL; child = child->next) {
    if (child->type != XML_ELEMENT_NODE) {
      continue;
    }
    if (current++ == index) {
      *out_child = elff_xml_node_to_handle(child);
      return ELFF_XML_OK;
    }
  }
  elff_xml_error_set(error,
                     elff_xml_node_line_number(parent_node),
                     0,
                     "The XML element-child index is out of range.");
  return ELFF_XML_NOT_FOUND;
#endif
}

elff_xml_status
elff_xml_node_name(elff_xml_node_t node,
                   const char** out_name,
                   elff_xml_error* error)
{
  elff_xml_error_clear(error);
  if (out_name != NULL) {
    *out_name = NULL;
  }
  if (node == NULL || out_name == NULL) {
    elff_xml_error_set(error, 0, 0, "Invalid XML node-name argument.");
    return ELFF_XML_INVALID_ARGUMENT;
  }
#ifndef ELFF_USE_LIBXML2
  elff_xml_error_set(error, 0, 0, "ELFF was built without libxml2 support.");
  return ELFF_XML_UNAVAILABLE;
#else
  const xmlNode* element = elff_xml_node_from_handle(node);
  if (element->type != XML_ELEMENT_NODE || element->ns != NULL) {
    elff_xml_error_set(error,
                       elff_xml_node_line_number(element),
                       0,
                       "The selected node is not an unqualified element.");
    return ELFF_XML_INVALID_VALUE;
  }
  *out_name = (const char*)element->name;
  return ELFF_XML_OK;
#endif
}

elff_xml_status
elff_xml_node_line(elff_xml_node_t node,
                   size_t* out_line,
                   elff_xml_error* error)
{
  elff_xml_error_clear(error);
  if (out_line != NULL) {
    *out_line = 0;
  }
  if (node == NULL || out_line == NULL) {
    elff_xml_error_set(error, 0, 0, "Invalid XML node-line argument.");
    return ELFF_XML_INVALID_ARGUMENT;
  }
#ifndef ELFF_USE_LIBXML2
  elff_xml_error_set(error, 0, 0, "ELFF was built without libxml2 support.");
  return ELFF_XML_UNAVAILABLE;
#else
  *out_line = elff_xml_node_line_number(elff_xml_node_from_handle(node));
  return ELFF_XML_OK;
#endif
}

elff_xml_status
elff_xml_node_attribute_string(elff_xml_node_t node,
                               const char* attribute_name,
                               char* value,
                               size_t value_capacity,
                               size_t* required_capacity,
                               elff_xml_error* error)
{
  elff_xml_error_clear(error);

  if (required_capacity != NULL) {
    *required_capacity = 0;
  }
  if (node == NULL || attribute_name == NULL || attribute_name[0] == '\0' ||
      (value == NULL && value_capacity != 0)) {
    elff_xml_error_set(error, 0, 0, "Invalid XML attribute argument.");
    return ELFF_XML_INVALID_ARGUMENT;
  }

#ifndef ELFF_USE_LIBXML2
  (void)value;
  (void)value_capacity;
  elff_xml_error_set(error, 0, 0, "ELFF was built without libxml2 support.");
  return ELFF_XML_UNAVAILABLE;
#else
  const xmlNode* element = elff_xml_node_from_handle(node);
  xmlChar* property = xmlGetProp((xmlNode*)element,
                                 (const xmlChar*)attribute_name);
  size_t needed;

  if (property == NULL) {
    char message[ELFF_XML_ERROR_MESSAGE_CAPACITY];
    (void)snprintf(message,
                   sizeof(message),
                   "Required attribute '%s' is missing.",
                   attribute_name);
    elff_xml_error_set(error, elff_xml_node_line_number(element), 0, message);
    return ELFF_XML_NOT_FOUND;
  }

  needed = strlen((const char*)property) + 1;
  if (required_capacity != NULL) {
    *required_capacity = needed;
  }

  if (value == NULL || value_capacity < needed) {
    if (value != NULL && value_capacity > 0) {
      value[0] = '\0';
    }
    xmlFree(property);
    elff_xml_error_set(error,
                       elff_xml_node_line_number(element),
                       0,
                       "The XML attribute output buffer is too small.");
    return ELFF_XML_BUFFER_TOO_SMALL;
  }

  (void)memcpy(value, property, needed);
  xmlFree(property);
  return ELFF_XML_OK;
#endif
}

elff_xml_status
elff_xml_node_attribute_int(elff_xml_node_t node,
                            const char* attribute_name,
                            int* value,
                            elff_xml_error* error)
{
  char buffer[128];
  elff_xml_status status;

  if (value == NULL) {
    elff_xml_error_clear(error);
    elff_xml_error_set(error, 0, 0, "Invalid integer attribute output.");
    return ELFF_XML_INVALID_ARGUMENT;
  }

  status = elff_xml_node_attribute_string(
    node, attribute_name, buffer, sizeof(buffer), NULL, error);
  if (status != ELFF_XML_OK) {
    return status;
  }

#ifndef ELFF_USE_LIBXML2
  return ELFF_XML_UNAVAILABLE;
#else
  char* end = NULL;
  const char* start;
  long parsed;

  start = elff_xml_skip_space(buffer);
  errno = 0;
  parsed = strtol(start, &end, 10);
  if (start == end || errno == ERANGE || parsed < INT_MIN || parsed > INT_MAX ||
      !elff_xml_only_trailing_space(end)) {
    elff_xml_error_set(error,
                       elff_xml_node_line_number(elff_xml_node_from_handle(node)),
                       0,
                       "The XML attribute is not a valid integer.");
    return ELFF_XML_INVALID_VALUE;
  }

  *value = (int)parsed;
  return ELFF_XML_OK;
#endif
}

elff_xml_status
elff_xml_node_attribute_double(elff_xml_node_t node,
                               const char* attribute_name,
                               double* value,
                               elff_xml_error* error)
{
  char buffer[128];
  elff_xml_status status;

  if (value == NULL) {
    elff_xml_error_clear(error);
    elff_xml_error_set(error, 0, 0, "Invalid floating-point attribute output.");
    return ELFF_XML_INVALID_ARGUMENT;
  }

  status = elff_xml_node_attribute_string(
    node, attribute_name, buffer, sizeof(buffer), NULL, error);
  if (status != ELFF_XML_OK) {
    return status;
  }

#ifndef ELFF_USE_LIBXML2
  return ELFF_XML_UNAVAILABLE;
#else
  char* end = NULL;
  const char* start;
  double parsed;

  start = elff_xml_skip_space(buffer);
  errno = 0;
  parsed = strtod(start, &end);
  if (start == end || errno == ERANGE || !isfinite(parsed) ||
      !elff_xml_only_trailing_space(end)) {
    elff_xml_error_set(error,
                       elff_xml_node_line_number(elff_xml_node_from_handle(node)),
                       0,
                       "The XML attribute is not a valid finite number.");
    return ELFF_XML_INVALID_VALUE;
  }

  *value = parsed;
  return ELFF_XML_OK;
#endif
}

const char*
elff_xml_status_string(elff_xml_status status)
{
  switch (status) {
    case ELFF_XML_OK:
      return "ok";
    case ELFF_XML_IO_ERROR:
      return "I/O error";
    case ELFF_XML_MALFORMED_DOCUMENT:
      return "malformed document";
    case ELFF_XML_INVALID_ARGUMENT:
      return "invalid argument";
    case ELFF_XML_NOT_FOUND:
      return "not found";
    case ELFF_XML_NOT_UNIQUE:
      return "not unique";
    case ELFF_XML_INVALID_VALUE:
      return "invalid value";
    case ELFF_XML_BUFFER_TOO_SMALL:
      return "buffer too small";
    case ELFF_XML_INTERNAL_ERROR:
      return "internal error";
    case ELFF_XML_UNAVAILABLE:
      return "libxml2 unavailable";
    default:
      return "unknown status";
  }
}
