#include "elff/io/C/xmlCase.h"

#include <math.h>

elff_xml_status
elff_xml_node_attribute_double_or(elff_xml_node_t node,
                                  const char* attribute_name,
                                  double default_value,
                                  double* value,
                                  elff_xml_error* error)
{
  elff_xml_status status;

  if (value == NULL || !isfinite(default_value)) {
    elff_xml_error_clear(error);
    elff_xml_error_set(error,
                       0,
                       0,
                       "Invalid optional floating-point attribute argument.");
    return ELFF_XML_INVALID_ARGUMENT;
  }

  status =
    elff_xml_node_attribute_double(node, attribute_name, value, error);
  if (status == ELFF_XML_NOT_FOUND) {
    *value = default_value;
    elff_xml_error_clear(error);
    return ELFF_XML_OK;
  }
  return status;
}

elff_xml_status
elff_xml_node_attribute_xyz(elff_xml_node_t node,
                            double default_z,
                            double value[3],
                            elff_xml_error* error)
{
  elff_xml_status status;

  if (value == NULL) {
    elff_xml_error_clear(error);
    elff_xml_error_set(error, 0, 0, "Invalid XYZ attribute output.");
    return ELFF_XML_INVALID_ARGUMENT;
  }

  status = elff_xml_node_attribute_double(node, "x", &value[0], error);
  if (status != ELFF_XML_OK) {
    return status;
  }
  status = elff_xml_node_attribute_double(node, "y", &value[1], error);
  if (status != ELFF_XML_OK) {
    return status;
  }
  return elff_xml_node_attribute_double_or(
    node, "z", default_z, &value[2], error);
}

static elff_xml_status
elff_xml_validation_status(elff_xml_dtd_status status)
{
  switch (status) {
    case ELFF_XML_DTD_VALID:
      return ELFF_XML_OK;
    case ELFF_XML_DTD_INVALID_ARGUMENT:
      return ELFF_XML_INVALID_ARGUMENT;
    case ELFF_XML_DTD_UNAVAILABLE:
      return ELFF_XML_UNAVAILABLE;
    case ELFF_XML_DTD_MALFORMED_XML:
      return ELFF_XML_MALFORMED_DOCUMENT;
    case ELFF_XML_DTD_SUBTREE_NOT_FOUND:
      return ELFF_XML_NOT_FOUND;
    case ELFF_XML_DTD_SUBTREE_NOT_UNIQUE:
      return ELFF_XML_NOT_UNIQUE;
    case ELFF_XML_DTD_INVALID_DOCUMENT:
      return ELFF_XML_INVALID_VALUE;
    case ELFF_XML_DTD_MALFORMED_DTD:
    case ELFF_XML_DTD_INTERNAL_ERROR:
    default:
      return ELFF_XML_INTERNAL_ERROR;
  }
}

elff_xml_status
elff_xml_case_load_section_file(const char* path,
                                const char* container_element_name,
                                const char* section_element_name,
                                const char* dtd,
                                size_t dtd_size,
                                elff_xml_document_t* out_document,
                                elff_xml_node_t* out_section,
                                elff_xml_error* error)
{
  elff_xml_document_t document = NULL;
  elff_xml_node_t root = NULL;
  elff_xml_node_t section = NULL;
  elff_xml_status status;
  elff_xml_dtd_status validation_status;
  char diagnostic[ELFF_XML_ERROR_MESSAGE_CAPACITY] = { 0 };
  size_t line = 0;

  elff_xml_error_clear(error);
  if (out_document != NULL) {
    *out_document = NULL;
  }
  if (out_section != NULL) {
    *out_section = NULL;
  }
  if (path == NULL || container_element_name == NULL ||
      container_element_name[0] == '\0' || section_element_name == NULL ||
      section_element_name[0] == '\0' || dtd == NULL ||
      out_document == NULL || out_section == NULL) {
    elff_xml_error_set(error, 0, 0, "Invalid XML case-section argument.");
    return ELFF_XML_INVALID_ARGUMENT;
  }

  status = elff_xml_document_load_file(path, &document, error);
  if (status != ELFF_XML_OK) {
    return status;
  }
  status = elff_xml_document_root(
    document, container_element_name, &root, error);
  if (status != ELFF_XML_OK) {
    goto cleanup;
  }
  status = elff_xml_node_unique_child(
    root, section_element_name, &section, error);
  if (status != ELFF_XML_OK) {
    goto cleanup;
  }

  validation_status = elff_xml_dtd_validate_node(
    section, dtd, dtd_size, diagnostic, sizeof(diagnostic));
  status = elff_xml_validation_status(validation_status);
  if (status != ELFF_XML_OK) {
    (void)elff_xml_node_line(section, &line, NULL);
    elff_xml_error_set(error,
                       line,
                       0,
                       diagnostic[0] != '\0'
                         ? diagnostic
                         : elff_xml_dtd_status_string(validation_status));
    goto cleanup;
  }

  *out_document = document;
  *out_section = section;
  return ELFF_XML_OK;

cleanup:
  elff_xml_document_destroy(document);
  return status;
}
