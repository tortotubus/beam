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

/** Render numeric ${expression} attributes and leaf text into an independent in-memory Case
 * document. Bare names and parameters.NAME query parameters. Definitions are retained unchanged for provenance.
 * V1 requires Parameters when using expressions; v2 permits constant expressions. Expressions
 * must occupy the whole value (apart from whitespace). Source is unchanged.
 * Legacy cases without definitions or expressions are copied unchanged.
 * On failure *out_document is null; callers own successful documents. */
elff_xml_status elff_xml_case_render(elff_xml_document_t source,
                                    elff_xml_document_t* out_document,
                                    elff_xml_error* error);

/** Serialize a document; this is never used to feed downstream readers. */
elff_xml_status elff_xml_document_write_file(elff_xml_document_t document,
                                             const char* path,
                                             elff_xml_error* error);

/** Load then render in memory. If write_diagnostic is nonzero, also write
 * NAME.rendered.xml beside NAME.xml (otherwise append .rendered.xml).
 * Diagnostic write failures warn on stderr but do not fail this operation.
 * MPI callers must enable diagnostics on only one rank. */
elff_xml_status elff_xml_case_load_rendered_file(const char* path,
                                                int write_diagnostic,
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

/** @brief Read finite numeric text from a leaf element; expressions must be resolved first. */
elff_xml_status
elff_xml_node_text_double(elff_xml_node_t node, double* value, elff_xml_error* error);

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

/** @brief Validate the complete Case using the embedded DTD selected by formatVersion.
 * Supports versions 1 and 2. Accepts source expressions or rendered values; numeric and contextual
 * rules remain the responsibility of the parameter and section readers.
 * Does not mutate the document or load external DTDs.
 */
elff_xml_status
elff_xml_case_validate(elff_xml_document_t document, elff_xml_error* error);

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
