#pragma once

#include <stddef.h>

#ifdef __cplusplus
extern "C"
{
#endif

#define ELFF_XML_ERROR_MESSAGE_CAPACITY 512

/** Opaque owner of one parsed, read-only XML document. */
typedef struct elff_xml_document_handle* elff_xml_document_t;

/** Borrowed element handle, valid while its document remains alive. */
typedef const struct elff_xml_node_handle* elff_xml_node_t;

typedef enum
{
  ELFF_XML_OK = 0,
  ELFF_XML_IO_ERROR,
  ELFF_XML_MALFORMED_DOCUMENT,
  ELFF_XML_INVALID_ARGUMENT,
  ELFF_XML_NOT_FOUND,
  ELFF_XML_NOT_UNIQUE,
  ELFF_XML_INVALID_VALUE,
  ELFF_XML_BUFFER_TOO_SMALL,
  ELFF_XML_INTERNAL_ERROR,
  ELFF_XML_UNAVAILABLE
} elff_xml_status;

typedef struct
{
  size_t line;
  size_t column;
  char message[ELFF_XML_ERROR_MESSAGE_CAPACITY];
} elff_xml_error;

/** Clear a diagnostic. Null is accepted. */
void
elff_xml_error_clear(elff_xml_error* error);

/** Set a diagnostic, trimming trailing newlines. Null is accepted. */
void
elff_xml_error_set(elff_xml_error* error,
                   size_t line,
                   size_t column,
                   const char* message);

elff_xml_status
elff_xml_document_load_memory(const char* xml,
                              size_t xml_size,
                              elff_xml_document_t* out_document,
                              elff_xml_error* error);

elff_xml_status
elff_xml_document_load_file(const char* path,
                            elff_xml_document_t* out_document,
                            elff_xml_error* error);

void
elff_xml_document_destroy(elff_xml_document_t document);

elff_xml_status
elff_xml_document_root(elff_xml_document_t document,
                       const char* expected_name,
                       elff_xml_node_t* out_root,
                       elff_xml_error* error);

elff_xml_status
elff_xml_node_unique_child(elff_xml_node_t parent,
                           const char* child_name,
                           elff_xml_node_t* out_child,
                           elff_xml_error* error);

elff_xml_status
elff_xml_node_element_child_count(elff_xml_node_t parent,
                                  size_t* out_count,
                                  elff_xml_error* error);

elff_xml_status
elff_xml_node_element_child_at(elff_xml_node_t parent,
                               size_t index,
                               elff_xml_node_t* out_child,
                               elff_xml_error* error);

elff_xml_status
elff_xml_node_name(elff_xml_node_t node,
                   const char** out_name,
                   elff_xml_error* error);

elff_xml_status
elff_xml_node_line(elff_xml_node_t node,
                   size_t* out_line,
                   elff_xml_error* error);

elff_xml_status
elff_xml_node_attribute_string(elff_xml_node_t node,
                               const char* attribute_name,
                               char* value,
                               size_t value_capacity,
                               size_t* required_capacity,
                               elff_xml_error* error);

elff_xml_status
elff_xml_node_attribute_int(elff_xml_node_t node,
                            const char* attribute_name,
                            int* value,
                            elff_xml_error* error);

elff_xml_status
elff_xml_node_attribute_double(elff_xml_node_t node,
                               const char* attribute_name,
                               double* value,
                               elff_xml_error* error);

/** Read a finite double, using default_value only when the attribute is absent. */
elff_xml_status
elff_xml_node_attribute_double_or(elff_xml_node_t node,
                                  const char* attribute_name,
                                  double default_value,
                                  double* value,
                                  elff_xml_error* error);

/** Read x and y, and read z or use default_z when z is absent. */
elff_xml_status
elff_xml_node_attribute_xyz(elff_xml_node_t node,
                            double default_z,
                            double value[3],
                            elff_xml_error* error);

const char*
elff_xml_status_string(elff_xml_status status);

typedef enum
{
  ELFF_XML_DTD_VALID = 0,
  ELFF_XML_DTD_INVALID_DOCUMENT,
  ELFF_XML_DTD_MALFORMED_XML,
  ELFF_XML_DTD_MALFORMED_DTD,
  ELFF_XML_DTD_INVALID_ARGUMENT,
  ELFF_XML_DTD_INTERNAL_ERROR,
  ELFF_XML_DTD_UNAVAILABLE,
  ELFF_XML_DTD_SUBTREE_NOT_FOUND,
  ELFF_XML_DTD_SUBTREE_NOT_UNIQUE
} elff_xml_dtd_status;

elff_xml_dtd_status
elff_xml_dtd_validate_memory(const char* xml,
                             size_t xml_size,
                             const char* dtd,
                             size_t dtd_size,
                             char* diagnostic,
                             size_t diagnostic_capacity);

/** Validate only the selected node; its ancestors and siblings are ignored. */
elff_xml_dtd_status
elff_xml_dtd_validate_node(elff_xml_node_t node,
                           const char* dtd,
                           size_t dtd_size,
                           char* diagnostic,
                           size_t diagnostic_capacity);

elff_xml_dtd_status
elff_xml_dtd_validate_subtree_memory(const char* xml,
                                     size_t xml_size,
                                     const char* container_element_name,
                                     const char* subtree_element_name,
                                     const char* dtd,
                                     size_t dtd_size,
                                     char* diagnostic,
                                     size_t diagnostic_capacity);

const char*
elff_xml_dtd_status_string(elff_xml_dtd_status status);

/**
 * Load a case file, select one unique top-level section, and validate only that
 * section against its owner's DTD. On success the document owns out_section;
 * the caller must destroy it with elff_xml_document_destroy().
 */
elff_xml_status
elff_xml_case_load_section_file(const char* path,
                                const char* container_element_name,
                                const char* section_element_name,
                                const char* dtd,
                                size_t dtd_size,
                                elff_xml_document_t* out_document,
                                elff_xml_node_t* out_section,
                                elff_xml_error* error);

#ifdef __cplusplus
}
#endif
