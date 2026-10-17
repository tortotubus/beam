#pragma once

#include <stddef.h>

#ifdef __cplusplus
extern "C"
{
#endif

#define ELFF_XML_ERROR_MESSAGE_CAPACITY 512

/** @brief Opaque owner of one parsed, read-only XML document. */
typedef struct elff_xml_document_handle* elff_xml_document_t;

/**
 * @brief Borrowed element handle.
 *
 * A node remains valid only while its owning @ref elff_xml_document_t exists.
 */
typedef const struct elff_xml_node_handle* elff_xml_node_t;

/** @brief Result of parsing or reading values from an XML document. */
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

/** @brief Optional diagnostic populated by XML document operations. */
typedef struct
{
  size_t line;
  size_t column;
  char message[ELFF_XML_ERROR_MESSAGE_CAPACITY];
} elff_xml_error;

/**
 * @brief Parse an XML document from memory.
 *
 * The byte input need not be null-terminated. External entities and network
 * access are disabled. On failure, @p out_document is set to null.
 */
elff_xml_status
elff_xml_document_load_memory(const char* xml,
                              size_t xml_size,
                              elff_xml_document_t* out_document,
                              elff_xml_error* error);

/**
 * @brief Parse an XML document from a local file.
 *
 * External entities and network access are disabled. On failure,
 * @p out_document is set to null.
 */
elff_xml_status
elff_xml_document_load_file(const char* path,
                            elff_xml_document_t* out_document,
                            elff_xml_error* error);

/** @brief Destroy a parsed document and invalidate all of its node handles. */
void
elff_xml_document_destroy(elff_xml_document_t document);

/**
 * @brief Get the document element and require an unqualified element name.
 */
elff_xml_status
elff_xml_document_root(elff_xml_document_t document,
                       const char* expected_name,
                       elff_xml_node_t* out_root,
                       elff_xml_error* error);

/**
 * @brief Select exactly one direct, unqualified child element by name.
 */
elff_xml_status
elff_xml_node_unique_child(elff_xml_node_t parent,
                           const char* child_name,
                           elff_xml_node_t* out_child,
                           elff_xml_error* error);

/**
 * @brief Copy an attribute value into a caller-owned buffer.
 *
 * @p required_capacity receives the number of bytes needed including the null
 * terminator. Passing a null buffer with zero capacity may be used to query
 * that size and returns @ref ELFF_XML_BUFFER_TOO_SMALL for a present value.
 */
elff_xml_status
elff_xml_node_attribute_string(elff_xml_node_t node,
                               const char* attribute_name,
                               char* value,
                               size_t value_capacity,
                               size_t* required_capacity,
                               elff_xml_error* error);

/** @brief Read a complete base-10 integer attribute value. */
elff_xml_status
elff_xml_node_attribute_int(elff_xml_node_t node,
                            const char* attribute_name,
                            int* value,
                            elff_xml_error* error);

/** @brief Read a complete, finite floating-point attribute value. */
elff_xml_status
elff_xml_node_attribute_double(elff_xml_node_t node,
                               const char* attribute_name,
                               double* value,
                               elff_xml_error* error);

/** @brief Return a stable, human-readable name for an XML status. */
const char*
elff_xml_status_string(elff_xml_status status);

#ifdef __cplusplus
}
#endif
