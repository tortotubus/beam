#pragma once

#include "xmlCase.h"

#ifdef __cplusplus
extern "C" {
#endif

/** Owned, resolved numeric parameters and nondimensional coefficients. */
typedef struct elff_xml_parameters_handle *elff_xml_parameters_t;

/**
 * Explicitly read Parameters and optional Nondimensionalization children of
 * Case. Parameters require name and exactly one of value/expression;
 * Coefficients require name and expression. Expressions reference parameters
 * only, with forward references allowed. Names use [A-Za-z_][A-Za-z0-9_]*; pi
 * is reserved. Grammar: decimal/scientific numbers, pi, sqrt(...), + - * / ^
 * and parentheses. Powers associate right and bind more tightly than unary
 * signs. All values are finite doubles; no units, substitution or model binding
 * occurs. Limits: 4096 bytes per attribute, 4096 entries per section, and 128
 * levels each of expression nesting and parameter dependency recursion. On
 * failure *out_context is null. On success the context owns all its data and
 * does not depend on the source document's lifetime. Uses elff_xml_error for
 * diagnostics, including the offending definition's XML line when available.
 */
elff_xml_status elff_xml_parameters_create(elff_xml_node_t case_node,
                                           elff_xml_parameters_t *out_context,
                                           elff_xml_error *error);

/** Query an exact name. Missing names return NOT_FOUND; output is unchanged on
 * failure. */
elff_xml_status elff_xml_parameters_get(elff_xml_parameters_t context,
                                        const char *name, double *value,
                                        elff_xml_error *error);
elff_xml_status
elff_xml_parameters_get_coefficient(elff_xml_parameters_t context,
                                    const char *name, double *value,
                                    elff_xml_error *error);

/** Evaluate a render expression: bare names and parameters.NAME refer to
 * parameters; coefficients.NAME refers to a coefficient. Output is unchanged
 * on failure. No XML is modified. */
elff_xml_status elff_xml_parameters_evaluate(elff_xml_parameters_t context,
                                            const char *expression, double *value,
                                            elff_xml_error *error);

/** Null is accepted. */
void elff_xml_parameters_destroy(elff_xml_parameters_t context);

#ifdef __cplusplus
}
#endif
