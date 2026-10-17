#pragma once

#include "elff/io/C/xmlCaseDocument.h"

#include <stddef.h>

#ifdef __cplusplus
extern "C"
{
#endif

/**
 * @brief Result of validating an XML document against a DTD.
 */
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

/**
 * @brief Validate an in-memory XML document against an in-memory DTD.
 *
 * Neither input needs to be null-terminated. External entities and network
 * access are disabled. When @p diagnostic is non-null and
 * @p diagnostic_capacity is nonzero, it receives a null-terminated error
 * description. A successful validation writes an empty diagnostic.
 *
 * @param xml XML document bytes
 * @param xml_size Number of bytes in @p xml
 * @param dtd DTD bytes
 * @param dtd_size Number of bytes in @p dtd
 * @param diagnostic Optional output buffer for an error description
 * @param diagnostic_capacity Size of @p diagnostic in bytes
 * @return Validation status
 */
elff_xml_dtd_status
elff_xml_dtd_validate_memory(const char* xml,
                             size_t xml_size,
                             const char* dtd,
                             size_t dtd_size,
                             char* diagnostic,
                             size_t diagnostic_capacity);

/**
 * @brief Validate a selected XML element as a standalone document.
 *
 * The node is copied from its owning parsed document before validation. The
 * handle remains borrowed and must belong to a live @ref elff_xml_document_t.
 * Siblings and ancestors are outside the scope of validation.
 *
 * @param node Borrowed document element to validate
 * @param dtd DTD bytes whose document element is @p node
 * @param dtd_size Number of bytes in @p dtd
 * @param diagnostic Optional output buffer for an error description
 * @param diagnostic_capacity Size of @p diagnostic in bytes
 * @return Validation status
 */
elff_xml_dtd_status
elff_xml_dtd_validate_node(elff_xml_node_t node,
                           const char* dtd,
                           size_t dtd_size,
                           char* diagnostic,
                           size_t diagnostic_capacity);

/**
 * @brief Validate one direct-child subtree against an in-memory DTD.
 *
 * The complete XML input is first parsed for well-formedness. Its unqualified
 * document element must be named @p container_element_name and must contain
 * exactly one direct, unqualified child named @p subtree_element_name. That
 * child is copied into a temporary standalone document and validated against
 * @p dtd. Sibling subtrees are therefore outside the scope of DTD validation.
 *
 * This intentionally implements a narrow selector rather than XPath. It is
 * suitable for documents such as `<Case><ELFF/> <Basilisk/></Case>`, where
 * each consumer owns and validates one section. Neither byte input needs to
 * be null-terminated; element names must be null-terminated. External
 * entities and network access are disabled.
 *
 * @param xml Complete XML document bytes
 * @param xml_size Number of bytes in @p xml
 * @param container_element_name Expected document-element name
 * @param subtree_element_name Direct-child element to select
 * @param dtd DTD bytes whose document element is the selected subtree
 * @param dtd_size Number of bytes in @p dtd
 * @param diagnostic Optional output buffer for an error description
 * @param diagnostic_capacity Size of @p diagnostic in bytes
 * @return Validation or subtree-selection status
 */
elff_xml_dtd_status
elff_xml_dtd_validate_subtree_memory(const char* xml,
                                     size_t xml_size,
                                     const char* container_element_name,
                                     const char* subtree_element_name,
                                     const char* dtd,
                                     size_t dtd_size,
                                     char* diagnostic,
                                     size_t diagnostic_capacity);

/**
 * @brief Return a stable, human-readable name for a validation status.
 *
 * @param status Validation status
 * @return Static string owned by the library
 */
const char*
elff_xml_dtd_status_string(elff_xml_dtd_status status);

#ifdef __cplusplus
}
#endif
