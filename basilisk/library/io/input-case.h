#pragma once

/* Header-only reader for the common <Case><Basilisk> settings.
 * No model construction, simulation globals, MPI calls or case dispatch.
 * read() borrows an already rendered document and copies its numeric/string
 * settings. It accepts unrelated caseType values and extension elements.
 * ImmersedBoundary is optional: its level defaults to Grid.maxLevel and
 * absent UniformForceDensity leaves the force disabled and zero.
 */
#include "elff/io/C/xmlCase.h"
#include <math.h>
#include <string.h>

typedef struct BasiliskCaseConfiguration
{
  double domain_size;
  double origin[3];
  double time_step;
  double end_time;
  double fluid_density;
  double dynamic_viscosity;
  double inlet_velocity[3];
  int min_level;
  int max_level;
  int immersed_boundary_level;
  int has_uniform_force_density;
  double uniform_force_density[3];
  char output_directory[256];
  double output_interval;
} BasiliskCaseConfiguration;

static elff_xml_status
basilisk_case_error (elff_xml_error * error, const char * message)
{
  elff_xml_error_set (error, 0, 0, message);
  return ELFF_XML_INVALID_VALUE;
}

static elff_xml_status
basilisk_case_read (elff_xml_document_t document,
                          BasiliskCaseConfiguration * configuration,
                          elff_xml_error * error)
{
  elff_xml_node_t basilisk = NULL, domain = NULL, time = NULL;
  elff_xml_node_t fluid = NULL, inlet = NULL, grid = NULL, immersed = NULL;
  elff_xml_node_t output = NULL, force = NULL;
  elff_xml_status status;

  if (configuration != NULL)
    memset (configuration, 0, sizeof (*configuration));
  if (document == NULL || configuration == NULL)
    return basilisk_case_error (error,
                                 "Invalid Basilisk case argument.");

  // Parse the Basilisk configuration without imposing a DTD. Required
  // children, typed attributes and numeric constraints are checked below.
  elff_xml_node_t root = NULL;
  status = elff_xml_document_root (document, "Case", &root, error);
  if (status == ELFF_XML_OK)
    status = elff_xml_node_unique_child (root, "Basilisk", &basilisk, error);
  if (status != ELFF_XML_OK)
    return status;

  if ((status = elff_xml_node_unique_child (
         basilisk, "Domain", &domain, error)) != ELFF_XML_OK ||
      (status = elff_xml_node_unique_child (
         basilisk, "Time", &time, error)) != ELFF_XML_OK ||
      (status = elff_xml_node_unique_child (
         basilisk, "Fluid", &fluid, error)) != ELFF_XML_OK ||
      (status = elff_xml_node_unique_child (
         fluid, "InletVelocity", &inlet, error)) != ELFF_XML_OK ||
      (status = elff_xml_node_unique_child (
         basilisk, "Grid", &grid, error)) != ELFF_XML_OK ||
      (status = elff_xml_node_unique_child (
         basilisk, "Output", &output, error)) != ELFF_XML_OK)
    return status;

  if ((status = elff_xml_node_attribute_double (
         domain, "size", &configuration->domain_size, error)) != ELFF_XML_OK ||
      (status = elff_xml_node_attribute_double (
         domain, "originX", &configuration->origin[0], error)) != ELFF_XML_OK ||
      (status = elff_xml_node_attribute_double (
         domain, "originY", &configuration->origin[1], error)) != ELFF_XML_OK ||
      (status = elff_xml_node_attribute_double_or (
         domain, "originZ", 0., &configuration->origin[2], error)) != ELFF_XML_OK ||
      (status = elff_xml_node_attribute_double (
         time, "step", &configuration->time_step, error)) != ELFF_XML_OK ||
      (status = elff_xml_node_attribute_double (
         time, "end", &configuration->end_time, error)) != ELFF_XML_OK ||
      (status = elff_xml_node_attribute_double (
         fluid, "density", &configuration->fluid_density, error)) != ELFF_XML_OK ||
      (status = elff_xml_node_attribute_double (
         fluid,
         "dynamicViscosity",
         &configuration->dynamic_viscosity,
         error)) != ELFF_XML_OK ||
      (status = elff_xml_node_attribute_xyz (
         inlet, 0., configuration->inlet_velocity, error)) != ELFF_XML_OK ||
      (status = elff_xml_node_attribute_int (
         grid, "minLevel", &configuration->min_level, error)) != ELFF_XML_OK ||
      (status = elff_xml_node_attribute_int (
         grid, "maxLevel", &configuration->max_level, error)) != ELFF_XML_OK ||
      (status = elff_xml_node_attribute_double (
         output, "interval", &configuration->output_interval, error)) != ELFF_XML_OK)
    return status;

  configuration->immersed_boundary_level = configuration->max_level;
  status = elff_xml_node_unique_child (
    basilisk, "ImmersedBoundary", &immersed, error);
  if (status == ELFF_XML_NOT_FOUND) {
    elff_xml_error_clear (error);
    immersed = NULL;
  }
  else if (status != ELFF_XML_OK)
    return status;
  if (immersed &&
      (status = elff_xml_node_attribute_int (
         immersed, "level", &configuration->immersed_boundary_level, error)) != ELFF_XML_OK)
    return status;

  status = elff_xml_node_attribute_string (
    output,
    "directory",
    configuration->output_directory,
    sizeof (configuration->output_directory),
    NULL,
    error);
  if (status != ELFF_XML_OK)
    return status;

  status = immersed ? elff_xml_node_unique_child (
    immersed, "UniformForceDensity", &force, error) : ELFF_XML_NOT_FOUND;
  if (status == ELFF_XML_NOT_FOUND) {
    configuration->has_uniform_force_density = 0;
    elff_xml_error_clear (error);
    status = ELFF_XML_OK;
  }
  else if (status == ELFF_XML_OK) {
    configuration->has_uniform_force_density = 1;
    status = elff_xml_node_attribute_xyz (
      force, 0., configuration->uniform_force_density, error);
    if (status != ELFF_XML_OK)
      return status;
  }
  else
    return status;

  if (!isfinite (configuration->domain_size) ||
      !isfinite (configuration->time_step) ||
      !isfinite (configuration->end_time) ||
      !isfinite (configuration->fluid_density) ||
      !isfinite (configuration->dynamic_viscosity) ||
      !isfinite (configuration->output_interval) ||
      !isfinite (configuration->origin[0]) ||
      !isfinite (configuration->origin[1]) ||
      !isfinite (configuration->origin[2]) ||
      !isfinite (configuration->inlet_velocity[0]) ||
      !isfinite (configuration->inlet_velocity[1]) ||
      !isfinite (configuration->inlet_velocity[2]) ||
      (configuration->has_uniform_force_density &&
       (!isfinite (configuration->uniform_force_density[0]) ||
        !isfinite (configuration->uniform_force_density[1]) ||
        !isfinite (configuration->uniform_force_density[2]))) ||
      !(configuration->domain_size > 0.) ||
      !(configuration->time_step > 0.) ||
      !(configuration->end_time > 0.) ||
      !(configuration->fluid_density > 0.) ||
      !(configuration->dynamic_viscosity >= 0.) ||
      configuration->min_level < 0 ||
      configuration->max_level < configuration->min_level ||
      configuration->immersed_boundary_level < 0 ||
      !(configuration->output_interval > 0.))
    return basilisk_case_error (error, "Invalid numeric value in <Basilisk> configuration.");

  return ELFF_XML_OK;
}

/** Load/render once and read the common settings. On success the caller owns
 * *document and may pass it to additional readers before destroying it.
 * On failure *document is null. write_diagnostic is explicit; MPI callers
 * should enable it on their designated writer only. No MPI is used here. */
static elff_xml_status
basilisk_case_load (const char *path, int write_diagnostic,
                   elff_xml_document_t *document,
                   BasiliskCaseConfiguration *configuration,
                   elff_xml_error *error)
{
  if (document) *document = NULL;
  if (!document || !configuration) {
    elff_xml_error_set (error, 0, 0, "Document output and configuration are required.");
    return ELFF_XML_INVALID_ARGUMENT;
  }
  elff_xml_status status = elff_xml_case_load_rendered_file (
    path, write_diagnostic, document, error);
  if (status == ELFF_XML_OK)
    status = basilisk_case_read (*document, configuration, error);
  if (status != ELFF_XML_OK) {
    elff_xml_document_destroy (*document);
    *document = NULL;
  }
  return status;
}
