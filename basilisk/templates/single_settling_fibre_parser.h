#pragma once
#include "library/io/input-case.h"

/* single_settling_fibre owns runtime initialization and its model requirements.
 * Include after library/elff/elff.h in this case's Basilisk translation unit.
 * Common settings are composed explicitly, not aliased to another case type.
 */
typedef struct SingleSettlingFibreCase {
  BasiliskCaseConfiguration basilisk;
} SingleSettlingFibreCase;

static int single_settling_fibre_case_init (const char *path, SingleSettlingFibreCase *configuration)
{
  elff_xml_document_t document = NULL;
  elff_xml_error error = {0};
  elff_xml_status status = configuration ? basilisk_case_load (
    path, pid() == 0, &document, &configuration->basilisk, &error) : ELFF_XML_INVALID_ARGUMENT;
  if (!configuration)
    elff_xml_error_set (&error, 0, 0, "Case configuration is required.");
  int success = status == ELFF_XML_OK;
#if _MPI && defined(ELFF_USE_MPI)
  int all_success = 0;
  MPI_Allreduce (&success, &all_success, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
#else
  int all_success = success;
#endif
  if (!all_success) {
    if (!success)
      fprintf (stderr, "single_settling_fibre XML error (line %zu): %s\n", error.line, error.message);
    elff_xml_document_destroy (document);
    return 1;
  }
  // Both readers consume this same document; runtime construction copies
  // what it needs, so no borrowed node survives document destruction.
  model_runtime_error runtime_error = {0};
  int rc = elff_init_document (document, &runtime_error);
  elff_xml_document_destroy (document);
  if (rc != 0) {
    fprintf (stderr, "ELFF model %zu (line %ld): %s\n",
             runtime_error.model_index, runtime_error.line, runtime_error.message);
    return rc;
  }
  ib_model_view_array_t models = elff_models ();
  if (models.count != 1) {
    fprintf (stderr, "single_settling_fibre requires exactly one fibre in <ELFF><Models>.\n");
    elff_deinit ();
    return 1;
  }
  return 0;
}
