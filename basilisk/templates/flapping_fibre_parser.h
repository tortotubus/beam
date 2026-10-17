#pragma once
#include "fibre_case_parser.h"

typedef FibreCaseConfiguration FlappingFibreConfiguration;
static elff_xml_status
flapping_fibre_case_load (const char *path, FlappingFibreConfiguration *config,
                          elff_xml_error *error)
{
  return fibre_case_load (path, config, error);
}
