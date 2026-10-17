#pragma once
#include "fibre_case_parser.h"

// Both consumers borrow this one rendered document. Diagnostic output is never
// read back and cannot affect the values supplied to either consumer.
static int fibre_case_init (const char *path, FibreCaseConfiguration *configuration) {
  elff_xml_document_t document = NULL;
  elff_xml_error error = {0};
  elff_xml_status status = elff_xml_case_load_rendered_file (
    path, pid() == 0, &document, &error);
  if (status == ELFF_XML_OK)
    status = fibre_case_read (document, configuration, &error);
  int success = status == ELFF_XML_OK;
#if _MPI && defined(ELFF_USE_MPI)
  int all_success = 0;
  MPI_Allreduce (&success, &all_success, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
#else
  int all_success = success;
#endif
  if (!all_success) {
    if (!success)
      fprintf (stderr, "Case XML error (line %zu): %s\n", error.line, error.message);
    elff_xml_document_destroy (document);
    return 1;
  }
  ib_runtime_error runtime_error = {0};
  int rc = elff_init_document (document, &runtime_error);
  elff_xml_document_destroy (document);
  if (rc != 0)
    fprintf (stderr, "ELFF model %zu (line %ld): %s\n",
             runtime_error.model_index, runtime_error.line, runtime_error.message);
  return rc;
}
