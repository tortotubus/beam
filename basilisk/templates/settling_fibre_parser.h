#pragma once

/** @brief Settling-specific controls layered on the common v2 configuration. */
typedef struct SettlingFibreCase {
  BasiliskXMLCase basilisk;
  int stats_interval;
} SettlingFibreCase;

/** @brief Default triply periodic Banaei settling configuration in simulation units. */
static int settling_fibre_defaults (BasiliskXMLCase* c) {
  basilisk_xml_case_defaults(c);
  for (int d = 0; d < 3; ++d) {
    c->grid_size[d] = 32.; c->grid_origin[d] = -16.; c->periodic_axis[d] = true;
  }
  c->grid_initial_level = c->grid_min_level = 4;
  c->grid_max_level = c->ib_level = 12;
  c->time_max_dt = 0.001; c->time_end = 10.;
  c->fluid_density = 1.; c->dynamic_viscosity = sqrt(0.1)/40.;
  c->force_density[1] = -0.1*pi/(4.*30.*30.);
  c->fluid_output.interval = c->solid_output.interval = 0.1;
  strcpy(c->output_directory, "settling_fibre_output");
  const BasiliskXMLField fields[] = {{"l2", false, 0}, {"omega_z", false, 0},
    {"p", false, 0}, {"u", true, 0}, {"ibmf", true, 0}};
  return basilisk_xml_default_fields(c, fields, 5, 0.001*sqrt(0.1));
}

/** @brief Validate reporting controls and the settling domain. */
static elff_xml_status settling_fibre_case_read (elff_xml_document_t document,
                                                SettlingFibreCase *configuration,
                                                elff_xml_error *error)
{
  double interval = 1.;
  elff_xml_status status = basilisk_xml_parameter (&configuration->basilisk,
                                                  "statsInterval", &interval, error);
  if (status != ELFF_XML_OK && status != ELFF_XML_NOT_FOUND) return status;
  elff_xml_error_clear (error);
  if (!isfinite (interval) || interval < 1. || interval > INT_MAX || floor (interval) != interval) {
    elff_xml_error_set (error, 0, 0, "Require a positive integer statsInterval fitting int.");
    return ELFF_XML_INVALID_VALUE;
  }
  configuration->stats_interval = (int) interval;
  for (int axis = 0; axis < 3; ++axis)
    if (!configuration->basilisk.periodic_axis[axis]) {
      elff_xml_error_set (error, 0, 0, "settling_fibre requires periodic boundaries on all three axes.");
      return ELFF_XML_INVALID_VALUE;
    }
  elff_xml_node_t root = NULL, elff = NULL;
  status = elff_xml_document_root (document, "Case", &root, error);
  if (status != ELFF_XML_OK) return status;
  if ((status = elff_xml_node_unique_child (root, "ELFF", &elff, error)) != ELFF_XML_OK)
    return status;
  int dimensions = 0;
  if ((status = elff_xml_node_attribute_int (root, "dimension", &dimensions, error)) != ELFF_XML_OK)
    return status;
  if (dimensions != 3) {
    elff_xml_error_set (error, 0, 0, "settling_fibre requires three-dimensional model geometry.");
    return ELFF_XML_INVALID_VALUE;
  }
  return ELFF_XML_OK;
}

/** @brief Render and read v2 settings before collective model construction. */
static int settling_fibre_case_init (const char *path, SettlingFibreCase *configuration)
{
  int defaults_ok = !settling_fibre_defaults(&configuration->basilisk);
  if (!generic_case_agree(defaults_ok)) {
    basilisk_xml_case_destroy(&configuration->basilisk); return 1;
  }
  elff_xml_document_t document = NULL;
  elff_xml_error error = {0};
  int success = elff_xml_case_load_rendered_file (path, pid() == 0, &document, &error) == ELFF_XML_OK;
  if (success)
    success = !basilisk_xml_case_read (document, "settling_fibre", &configuration->basilisk, &error);
  if (success)
    success = settling_fibre_case_read (document, configuration, &error) == ELFF_XML_OK;
  if (!generic_case_agree (success)) {
    if (!success) fprintf (stderr, "settling_fibre XML error (line %zu): %s\n", error.line, error.message);
    elff_xml_document_destroy (document);
    basilisk_xml_case_destroy (&configuration->basilisk);
    return 1;
  }
  model_runtime_error runtime_error = {0};
  int status = elff_init_document (document, &runtime_error);
  elff_xml_document_destroy (document);
  if (status) {
    fprintf (stderr, "ELFF model %zu (line %ld): %s\n", runtime_error.model_index,
             runtime_error.line, runtime_error.message);
    basilisk_xml_case_destroy (&configuration->basilisk);
    return status;
  }
  ib_model_view_array_t models = elff_models ();
  if (!models.count || models.count > INT_MAX/7) {
    fprintf (stderr, "settling_fibre requires a nonempty population fitting mesh/statistics indices.\n");
    elff_deinit ();
    basilisk_xml_case_destroy (&configuration->basilisk);
    return 1;
  }
  return 0;
}
