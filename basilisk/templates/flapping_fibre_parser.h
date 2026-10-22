#pragma once

/** @brief Flapping-specific initialization layered on the common v2 reader. */
typedef struct {
  BasiliskXMLCase basilisk;
  double initial_velocity;
} FlappingFibreCase;

/** @brief Read v2 configuration, agree on parsing success, and construct models.
 * inletSpeed is an optional resolved coefficient; its default is one.
 * Shared XML settings and named values remain owned until case cleanup.
 */
static int flapping_fibre_case_init (const char* path, FlappingFibreCase* configuration) {
  basilisk_xml_case_defaults(&configuration->basilisk);
  configuration->initial_velocity = 1.;
  elff_xml_document_t document = NULL;
  elff_xml_error error = {0};
  int success = elff_xml_case_load_rendered_file(path, pid() == 0, &document, &error) == ELFF_XML_OK;
  if (success) success = !basilisk_xml_case_read(document, "flapping_fibre",
                                                &configuration->basilisk, &error);
  if (success) {
    elff_xml_status status = basilisk_xml_coefficient(&configuration->basilisk,
      "inletSpeed", &configuration->initial_velocity, &error);
    success = status == ELFF_XML_OK || status == ELFF_XML_NOT_FOUND;
    if (status == ELFF_XML_NOT_FOUND) elff_xml_error_clear(&error);
  }
  if (!generic_case_agree(success)) {
    if (!success) fprintf(stderr, "flapping_fibre XML error (line %zu): %s\n", error.line, error.message);
    elff_xml_document_destroy(document);
    basilisk_xml_case_destroy(&configuration->basilisk);
    return 1;
  }
  model_runtime_error runtime_error = {0};
  int status = elff_init_document(document, &runtime_error);
  elff_xml_document_destroy(document);
  if (status) {
    fprintf(stderr, "ELFF model %zu (line %ld): %s\n",
            runtime_error.model_index, runtime_error.line, runtime_error.message);
    basilisk_xml_case_destroy(&configuration->basilisk);
    return status;
  }
  ib_model_view_array_t models = elff_models();
  if (!models.count || models.count > INT_MAX) {
    fprintf(stderr, "flapping_fibre requires a nonempty model population fitting mesh indices.\n");
    elff_deinit();
    basilisk_xml_case_destroy(&configuration->basilisk);
    return 1;
  }
  return 0;
}
