#pragma once

#include "library/ibm/IBMeshManager.h"
#include "library/elff/elff.h"
#include "library/ibm/IBOutput.h"
#include "library/io/output-vtk.h"
#include "elff/io/C/xmlCase.h"
#include "elff/io/C/xmlCaseParameters.h"

#include <math.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/**
 * Type Definitions
 */

/** @brief Boundary types supported by the generic case configuration.
 * Periodicity is configured per grid axis, not per field.
 */
typedef enum { BASILISK_BC_DIRICHLET, BASILISK_BC_NEUMANN } BasiliskBCType;

/** @brief Physical boundary sides, independent of Basilisk's internal numbering. */
typedef enum {
  BASILISK_BC_LEFT,
  BASILISK_BC_RIGHT,
  BASILISK_BC_BOTTOM,
  BASILISK_BC_TOP,
  BASILISK_BC_BACK,
  BASILISK_BC_FRONT
} BasiliskBCDirection;

/** @brief Owned constant boundary data; field_name is a scalar or Cartesian component.
 * Examples: p, pf, u.x, u.y. Whole vectors and face fields are not accepted.
 * source_line is retained for XML diagnostics. No XML nodes are borrowed.
 */
typedef struct {
  char field_name[128];
  double value;
  BasiliskBCType type;
  BasiliskBCDirection direction;
  size_t source_line;
} BasiliskBC;

/** @brief Simulation-time output cadence. Disabled schedules create no time events. */
typedef struct {
  bool enabled;
  double start, interval, end;
} BasiliskOutputSchedule;

/** @brief An owned output field name, bound after solver initialization. */
typedef struct {
  char name[128];
  bool vector_field;
  size_t source_line;
} BasiliskXMLField;

/** @brief Generic case settings, populated before run().
 * Boundary storage and named parameters/coefficients are owned; field lists
 * are borrowed unless owns_field_lists is set by the XML reader. Copies borrow the owned
 * storage and must not be destroyed separately. XML parsing is separate from field
 * binding and from the simulation events. Use basilisk_xml_case_defaults() before
 * populating this structure.
 */
typedef struct {
  char case_name[128];
  void (*initialize_fluid) (void);
  void (*prepare_fluid_output) (void);
  ib_kernel_1d kernel;
  double checkpoint_interval, checkpoint_end;
  bool owns_field_lists;
  BasiliskXMLField *fluid_fields, *solid_fields;
  size_t fluid_field_count, solid_field_count;
  elff_xml_parameters_t named_values;
  double grid_size[3], grid_origin[3];
  int grid_initial_level, grid_min_level, grid_max_level, ib_level;
  bool periodic_axis[3];
  double time_max_dt, time_end;
  double fluid_density, dynamic_viscosity;
  BasiliskBC* bcs;
  size_t bc_count;
  char output_directory[256];
  bool checkpoint_enable, restore_enable;
  BasiliskOutputSchedule fluid_output, solid_output;
  scalar* fluid_scalars;
  vector* fluid_vectors;
  IBscalar* solid_scalars;
  IBvector* solid_vectors;
  scalar* adapt_fields;
  double* adapt_tolerances;
  size_t adapt_tolerance_count;
  bool adapt_enable;
  double force_density[3];
  double force_relaxation;
  int force_iterations;
} BasiliskXMLCase;

/**
 * Function Declarations
 */

/** @brief Initialize defaults before reading or assigning configuration. */
static void basilisk_xml_case_defaults (BasiliskXMLCase* c);

/** @brief Release owned configuration storage after the simulation has finished. */
static void basilisk_xml_case_destroy (BasiliskXMLCase* c);

/** @brief Resolve and own Parameters and Nondimensionalization from a Case document.
 * Replaces named values only on success; the document may then be destroyed.
 * Case-specific readers choose recognized names, defaults and range checks.
 * Definitions do not automatically change simulation settings.
 * Both sections may be absent. Parameter expressions use the existing ELFF
 * evaluator; coefficients occupy a separate namespace. Call this reader before
 * the lookup helpers, even when no definitions are expected.
 * @code
 * double stats_interval = 10.;
 * elff_xml_status status = basilisk_xml_parameter(c, "statsInterval",
 *                                                &stats_interval, error);
 * // Keep the default on NOT_FOUND; report other failures, then validate that
 * // stats_interval is a positive integer before converting it.
 * @endcode
 */
static elff_xml_status basilisk_xml_read_named_values (elff_xml_document_t document,
                                                       BasiliskXMLCase* c,
                                                       elff_xml_error* error);

/** @brief Query a parameter by exact name; missing names return NOT_FOUND.
 * Output is unchanged on failure. A reader may supply its default on NOT_FOUND.
 */
static elff_xml_status basilisk_xml_parameter (const BasiliskXMLCase* c,
                                               const char* name,
                                               double* value,
                                               elff_xml_error* error);

/** @brief Query a resolved coefficient, independently of the parameter namespace.
 * Output is unchanged on failure; missing names return NOT_FOUND.
 */
static elff_xml_status basilisk_xml_coefficient (const BasiliskXMLCase* c,
                                                 const char* name,
                                                 double* value,
                                                 elff_xml_error* error);

/** @brief ... */
static int basilisk_xml_error (elff_xml_error* error, size_t line, const char* message);

/** @brief Check v2 grid selection and optional ELFF dimensions against this executable.
 * Call before constructing models or configuring the grid. Grid.type also fixes
 * the dimension for fluid-only cases. No simulation settings are modified.
 * Returns zero on success, otherwise a source-located diagnostic.
 */
static int basilisk_xml_validate_backend (elff_xml_document_t document,
                                          elff_xml_error* error);

/** @brief ... */
static int basilisk_xml_add_bc (BasiliskXMLCase* c,
                                const char* field,
                                BasiliskBCDirection side,
                                BasiliskBCType type,
                                double value,
                                size_t line,
                                elff_xml_error* error);
/** @brief ... */
static int basilisk_xml_validate_boundaries (const BasiliskXMLCase* c,
                                             elff_xml_error* error);

/** @brief ... */
static int basilisk_xml_bind_boundaries (const BasiliskXMLCase* c, elff_xml_error* error);

/** @brief Read v2 FluidBoundaryCondition children into a temporary owned list.
 * Numeric Value text is currently required; text expressions are not rendered
 * by the v1 renderer. Periodic declarations use field="*" and affect the axis.
 * This operation replaces the configured boundaries only after complete success.
 */
static int basilisk_xml_read_boundaries (elff_xml_node_t fluid,
                                         BasiliskXMLCase* c,
                                         elff_xml_error* error);

/** @brief Agree on setup success across ranks before entering later collectives. */
static int generic_case_agree (int success);

/** @brief Validate an enabled simulation-time output schedule. */
static bool generic_schedule_valid (BasiliskOutputSchedule schedule);

static int generic_ibm_case_run (const BasiliskXMLCase* configuration);

/** @brief Read a rendered v2 document; configuration must have defaults initialized.
 * On failure the previous configuration is unchanged. Parser-created field lists
 * are owned. Field names are bound after IB defaults have allocated their fields.
 */
static int basilisk_xml_case_read (elff_xml_document_t document,
                                   const char* expected_case,
                                   BasiliskXMLCase* configuration,
                                   elff_xml_error* error);
static int basilisk_xml_optional (elff_xml_node_t parent,
                                  const char* name,
                                  elff_xml_node_t* node,
                                  elff_xml_error* error);
static int
basilisk_xml_enabled (elff_xml_node_t node, bool* enabled, elff_xml_error* error);
static int basilisk_xml_schedule (elff_xml_node_t node,
                                  double end,
                                  BasiliskOutputSchedule* schedule,
                                  elff_xml_error* error);
static int basilisk_xml_output_fields (elff_xml_node_t node,
                                       bool solid,
                                       BasiliskXMLField** fields,
                                       size_t* count,
                                       elff_xml_error* error);
static int basilisk_xml_bind_output_fields (BasiliskXMLCase* c, elff_xml_error* error);

/**
 * Globals
 */

static int generic_case_status;
static double generic_end_time = 1.;
static double generic_fluid_start = 0., generic_fluid_interval = 1.,
              generic_fluid_end = 1.;
static double generic_solid_start = 0., generic_solid_interval = 1.,
              generic_solid_end = 1.;

static BasiliskXMLCase generic_case;

face vector generic_mu[], generic_alpha[];
scalar generic_rho[];

/** Macros */

#define BASILISK_XML_BIND(side, index)                                                   \
  case index:                                                                            \
    if (bc->type == BASILISK_BC_DIRICHLET) {                                             \
      s[side] = dirichlet (s.xml_bc_value[index]);                                       \
    } else {                                                                             \
      s[side] = neumann (s.xml_bc_value[index]);                                         \
    }                                                                                    \
    break;

/**
 * Events
 */

event properties (i++) {
  foreach_face () {
    generic_mu.x[] = fm.x[] * generic_case.dynamic_viscosity;
    generic_alpha.x[] = fm.x[] / generic_case.fluid_density;
  }
  foreach ()
    generic_rho[] = cm[] * generic_case.fluid_density;
}

event init (i = 0) {
  elff_xml_error error = {0};
  int success = !basilisk_xml_bind_boundaries (&generic_case, &error);
  if (!generic_case_agree (success)) {
    if (!success)
      fprintf (stderr,
               "Boundary configuration error (line %zu): %s\n",
               error.line,
               error.message);
    generic_case_status = 1;
    return 1;
  }
  ib_model_view_array_t models = elff_models ();
  for (size_t m = 0; m < models.count; ++m) {
    int mesh_id = ibmeshmanager_add_mesh ();
    success =
      elff_register_model (mesh_id, (int) m, models.data[m], generic_case.ib_level) == 0;
    if (!generic_case_agree (success)) {
      generic_case_status = 1;
      return 1;
    }
  }
  success = !basilisk_xml_bind_output_fields (&generic_case, &error);
  if (!generic_case_agree (success)) {
    if (!success)
      fprintf (stderr, "Output fields: %s\n", error.message);
    generic_case_status = 1;
    return 1;
  }
  bool restored =
    generic_case.restore_enable && restore_handler (generic_case.output_directory);
  if (!restored) {
    foreach ()
      foreach_dimension () u.x[] = 0.;
    if (generic_case.initialize_fluid)
      generic_case.initialize_fluid ();
#if TREE
    if (models.count)
      adapt_wavelet_ibm (NULL, NULL, 0, 1, all, true);
#endif
  }
  // Face material fields are derived and are not restored by the fluid dump.
  // Pressure boundary callbacks can divide by alpha, so rebuild them first.
  if (restored)
    event ("properties");
  boundary (all);
}

/* Inherited event executes before the solver advances structural models. */
event advance_interface (i++) {
  foreach_ibnode () {
    ibval (gravity.x) = generic_case.force_density[0];
    ibval (gravity.y) = generic_case.force_density[1];
#if dimension > 2
    ibval (gravity.z) = generic_case.force_density[2];
#endif
  }
}

event case_output_solids (t = generic_solid_start; t += generic_solid_interval;
                          t <= generic_solid_end) {
  if (generic_case.solid_output.enabled)
    output_hdf_pd (generic_case.solid_scalars,
                   generic_case.solid_vectors,
                   generic_case.output_directory);
}

event case_output_fluid (t = generic_fluid_start; t += generic_fluid_interval;
                         t <= generic_fluid_end) {
  if (generic_case.fluid_output.enabled) {
    if (generic_case.prepare_fluid_output)
      generic_case.prepare_fluid_output ();
#if TREE
    output_hdf_htg (generic_case.fluid_scalars,
                    generic_case.fluid_vectors,
                    generic_case.output_directory);
#else
    output_hdf_imagedata (generic_case.fluid_scalars,
                          generic_case.fluid_vectors,
                          generic_case.output_directory);
#endif
  }
}

#if TREE
event adapt (i++, last) {
  if (generic_case.adapt_enable)
    adapt_wavelet_ibm (generic_case.adapt_fields,
                       generic_case.adapt_tolerances,
                       generic_case.grid_max_level,
                       generic_case.grid_min_level);
}
#endif

event case_checkpoint (i++, last) {
  if (generic_case.checkpoint_enable) {
    checkpoint_configuration.checkpoint_on_sim_time =
      generic_case.checkpoint_interval > 0 && t <= generic_case.checkpoint_end;
    return checkpoint_handler (t, i, generic_case.output_directory);
  }
}

event case_stop (t = generic_end_time) {
  return 1;
}

event case_cleanup (i = end, last) {
  if (generic_case.owns_field_lists) {
    free (generic_case.fluid_scalars);
    free (generic_case.fluid_vectors);
    free (generic_case.solid_scalars);
    free (generic_case.solid_vectors);
    generic_case.fluid_scalars = NULL;
    generic_case.fluid_vectors = NULL;
    generic_case.solid_scalars = NULL;
    generic_case.solid_vectors = NULL;
  }
  ibmeshmanager_free ();
  elff_deinit ();
}

/** Function Definitions */

static int generic_case_agree (int success) {
#if _MPI
  int result = 0;
  MPI_Allreduce (&success, &result, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  return result;
#else
  return success;
#endif
}

static bool generic_schedule_valid (BasiliskOutputSchedule schedule) {
  return !schedule.enabled ||
         (isfinite (schedule.start) && isfinite (schedule.interval) &&
          isfinite (schedule.end) && schedule.start >= 0 && schedule.interval > 0 &&
          schedule.end >= 0);
}

static void basilisk_xml_case_defaults (BasiliskXMLCase* c) {
  memset (c, 0, sizeof (*c));
  for (int d = 0; d < 3; ++d)
    c->grid_size[d] = 1.;
  c->grid_initial_level = c->grid_min_level = 4;
  c->grid_max_level = c->ib_level = 6;
  c->time_max_dt = 0.001;
  c->time_end = 1.;
  c->fluid_density = 1.;
  c->force_relaxation = 0.8;
  c->kernel = peskin_cosine_kernel_1d;
  c->checkpoint_end = HUGE;
  c->force_iterations = 2;
  c->fluid_output.interval = c->solid_output.interval = 1.;
  c->fluid_output.end = c->solid_output.end = 1.;
  strcpy (c->output_directory, "output");
}

static void basilisk_xml_case_destroy (BasiliskXMLCase* c) {
  if (c->owns_field_lists) {
    free (c->fluid_scalars);
    free (c->fluid_vectors);
    free (c->solid_scalars);
    free (c->solid_vectors);
    free (c->adapt_fields);
    free (c->adapt_tolerances);
    c->fluid_scalars = NULL;
    c->fluid_vectors = NULL;
    c->solid_scalars = NULL;
    c->solid_vectors = NULL;
    c->adapt_fields = NULL;
    c->adapt_tolerances = NULL;
    free (c->fluid_fields);
    free (c->solid_fields);
    c->fluid_fields = c->solid_fields = NULL;
    c->fluid_field_count = c->solid_field_count = c->adapt_tolerance_count = 0;
    c->owns_field_lists = false;
  }
  elff_xml_parameters_destroy (c->named_values);
  c->named_values = NULL;
  free (c->bcs);
  c->bcs = NULL;
  c->bc_count = 0;
}

static elff_xml_status basilisk_xml_read_named_values (elff_xml_document_t document,
                                                       BasiliskXMLCase* c,
                                                       elff_xml_error* error) {
  if (!c) {
    elff_xml_error_set (error, 0, 0, "Case configuration is required.");
    return ELFF_XML_INVALID_ARGUMENT;
  }
  elff_xml_node_t root = NULL;
  elff_xml_status status = elff_xml_document_root (document, "Case", &root, error);
  if (status != ELFF_XML_OK)
    return status;
  elff_xml_parameters_t values = NULL;
  status = elff_xml_parameters_create (root, &values, error);
  if (status != ELFF_XML_OK)
    return status;
  elff_xml_parameters_destroy (c->named_values);
  c->named_values = values;
  return ELFF_XML_OK;
}

static elff_xml_status basilisk_xml_parameter (const BasiliskXMLCase* c,
                                               const char* name,
                                               double* value,
                                               elff_xml_error* error) {
  return elff_xml_parameters_get (c ? c->named_values : NULL, name, value, error);
}

static elff_xml_status basilisk_xml_coefficient (const BasiliskXMLCase* c,
                                                 const char* name,
                                                 double* value,
                                                 elff_xml_error* error) {
  return elff_xml_parameters_get_coefficient (
    c ? c->named_values : NULL, name, value, error);
}

static int basilisk_xml_error (elff_xml_error* error, size_t line, const char* message) {
  elff_xml_error_set (error, line, 0, message);
  return 1;
}

static int basilisk_xml_validate_backend (elff_xml_document_t document,
                                          elff_xml_error* error) {
  elff_xml_node_t root = NULL, basilisk = NULL, grid = NULL, elff = NULL;
  char type[64];
  size_t line = 0;
  if (elff_xml_document_root (document, "Case", &root, error) != ELFF_XML_OK ||
      elff_xml_node_unique_child (root, "Basilisk", &basilisk, error) != ELFF_XML_OK ||
      elff_xml_node_unique_child (basilisk, "Grid", &grid, error) != ELFF_XML_OK ||
      elff_xml_node_attribute_string (grid, "type", type, sizeof (type), NULL, error) !=
        ELFF_XML_OK)
    return 1;
  elff_xml_node_line (grid, &line, NULL);
  const char* expected = NULL;
  if (!strcmp (GRIDNAME, "Quadtree") && dimension == 2)
    expected = "quadtree";
  else if (!strcmp (GRIDNAME, "Octree") && dimension == 3)
    expected = "octree";
  else if (!strcmp (GRIDNAME, "Multigrid")) {
    if (dimension == 2)
      expected = "multigrid";
    else if (dimension == 3)
      expected = "multigrid3D";
  }
  if (!expected)
    return basilisk_xml_error (
      error, line, "Compiled grid backend is unsupported by the v2 XML schema.");
  if (strcmp (type, expected)) {
    char message[ELFF_XML_ERROR_MESSAGE_CAPACITY];
    snprintf (
      message,
      sizeof (message),
      "Grid type '%s' does not match this executable: expected '%s' (dimension %d).",
      type,
      expected,
      dimension);
    return basilisk_xml_error (error, line, message);
  }
  elff_xml_status status = elff_xml_node_unique_child (root, "ELFF", &elff, error);
  if (status == ELFF_XML_NOT_FOUND) {
    elff_xml_error_clear (error);
    return 0;
  }
  if (status != ELFF_XML_OK)
    return 1;
  int dimensions = 0;
  if (elff_xml_node_attribute_int (elff, "dimensions", &dimensions, error) != ELFF_XML_OK)
    return 1;
  if (dimensions != dimension) {
    char message[ELFF_XML_ERROR_MESSAGE_CAPACITY];
    elff_xml_node_line (elff, &line, NULL);
    snprintf (message,
              sizeof (message),
              "ELFF dimensions=%d does not match this executable's dimension %d.",
              dimensions,
              dimension);
    return basilisk_xml_error (error, line, message);
  }
  elff_xml_error_clear (error);
  return 0;
}

static int basilisk_xml_add_bc (BasiliskXMLCase* c,
                                const char* field,
                                BasiliskBCDirection side,
                                BasiliskBCType type,
                                double value,
                                size_t line,
                                elff_xml_error* error) {
  if (!field || !*field || strlen (field) >= sizeof (((BasiliskBC*) 0)->field_name) ||
      side < BASILISK_BC_LEFT || side > BASILISK_BC_FRONT ||
      (type != BASILISK_BC_DIRICHLET && type != BASILISK_BC_NEUMANN) || !isfinite (value))
    return basilisk_xml_error (
      error, line, "Invalid boundary name, side, type or value.");
  for (size_t j = 0; j < c->bc_count; ++j)
    if (c->bcs[j].direction == side && !strcmp (field, c->bcs[j].field_name))
      return basilisk_xml_error (error, line, "Duplicate field/side boundary condition.");
  if (c->bc_count >= (size_t) -1 / sizeof (BasiliskBC) - 1)
    return basilisk_xml_error (error, line, "Too many boundary conditions.");
  BasiliskBC* next = realloc (c->bcs, (c->bc_count + 1) * sizeof (BasiliskBC));
  if (!next)
    return basilisk_xml_error (error, line, "Cannot allocate boundary settings.");
  c->bcs = next;
  BasiliskBC* bc = &c->bcs[c->bc_count++];
  strcpy (bc->field_name, field);
  bc->direction = side;
  bc->type = type;
  bc->value = value;
  bc->source_line = line;
  return 0;
}

/* Values belong to the scalar, rather than a temporary XML node or loop variable.
 * qcc replaces s with the callback's scalar and generates homogeneous callbacks.
 */
attribute {
  double xml_bc_value[6];
}

static int basilisk_xml_validate_boundaries (const BasiliskXMLCase* c,
                                             elff_xml_error* error) {
  if (c->bc_count && !c->bcs)
    return basilisk_xml_error (error, 0, "Missing boundary storage.");
  for (int axis = dimension; axis < 3; ++axis)
    if (c->periodic_axis[axis])
      return basilisk_xml_error (
        error, 0, "Periodic axis is outside the compiled dimension.");
  for (size_t j = 0; j < c->bc_count; ++j) {
    const BasiliskBC* bc = &c->bcs[j];
    int axis = (int) bc->direction / 2;
    if (bc->direction < BASILISK_BC_LEFT || bc->direction > BASILISK_BC_FRONT ||
        axis >= dimension || c->periodic_axis[axis])
      return basilisk_xml_error (
        error, bc->source_line, "Boundary side is unavailable or its axis is periodic.");
    scalar s = lookup_field (bc->field_name);
    if (s.i < 0 || is_constant (s) || s.face)
      return basilisk_xml_error (
        error,
        bc->source_line,
        "Boundary requires an existing cell-centered scalar or vector component.");
    if (!isfinite (bc->value) ||
        (bc->type != BASILISK_BC_DIRICHLET && bc->type != BASILISK_BC_NEUMANN))
      return basilisk_xml_error (
        error, bc->source_line, "Invalid boundary value or type.");
    for (size_t k = 0; k < j; ++k)
      if (c->bcs[k].direction == bc->direction &&
          !strcmp (c->bcs[k].field_name, bc->field_name))
        return basilisk_xml_error (
          error, bc->source_line, "Duplicate field/side boundary.");
  }
  return 0;
}

static int basilisk_xml_bind_boundaries (const BasiliskXMLCase* c,
                                         elff_xml_error* error) {
  if (basilisk_xml_validate_boundaries (c, error))
    return 1;
  if (c->periodic_axis[0])
    periodic (right);
#if dimension > 1
  if (c->periodic_axis[1])
    periodic (top);
#endif
#if dimension > 2
  if (c->periodic_axis[2])
    periodic (front);
#endif
  for (size_t j = 0; j < c->bc_count; ++j) {
    const BasiliskBC* bc = &c->bcs[j];
    scalar s = lookup_field (bc->field_name);
    s.xml_bc_value[bc->direction] = bc->value;
    switch (bc->direction) {
      BASILISK_XML_BIND (left, BASILISK_BC_LEFT)
      BASILISK_XML_BIND (right, BASILISK_BC_RIGHT)
#if dimension > 1
      BASILISK_XML_BIND (bottom, BASILISK_BC_BOTTOM)
      BASILISK_XML_BIND (top, BASILISK_BC_TOP)
#endif
#if dimension > 2
      BASILISK_XML_BIND (back, BASILISK_BC_BACK)
      BASILISK_XML_BIND (front, BASILISK_BC_FRONT)
#endif
    }
#undef BASILISK_XML_BIND
  }
  return 0;
}

/** @brief Read v2 FluidBoundaryCondition children into a temporary owned list.
 * Numeric Value text is currently required; text expressions are not rendered
 * by the v1 renderer. Periodic declarations use field="*" and affect the axis.
 * This operation replaces the configured boundaries only after complete success.
 */
static int basilisk_xml_read_boundaries (elff_xml_node_t fluid,
                                         BasiliskXMLCase* c,
                                         elff_xml_error* error) {
  BasiliskXMLCase temporary;
  basilisk_xml_case_defaults (&temporary);
  size_t count = 0;
  if (elff_xml_node_element_child_count (fluid, &count, error) != ELFF_XML_OK)
    return 1;
  for (size_t j = 0; j < count; ++j) {
    elff_xml_node_t node = NULL, value_node = NULL;
    const char* name = NULL;
    char field[128], type[32], side[32];
    size_t line = 0;
    if (elff_xml_node_element_child_at (fluid, j, &node, error) != ELFF_XML_OK ||
        elff_xml_node_name (node, &name, error) != ELFF_XML_OK) {
      basilisk_xml_case_destroy (&temporary);
      return 1;
    }
    if (strcmp (name, "FluidBoundaryCondition"))
      continue;
    elff_xml_node_line (node, &line, NULL);
    if (elff_xml_node_attribute_string (
          node, "field", field, sizeof (field), NULL, error) != ELFF_XML_OK ||
        elff_xml_node_attribute_string (node, "type", type, sizeof (type), NULL, error) !=
          ELFF_XML_OK ||
        elff_xml_node_attribute_string (
          node, "direction", side, sizeof (side), NULL, error) != ELFF_XML_OK) {
      basilisk_xml_case_destroy (&temporary);
      return 1;
    }
    const char* sides[] = {"left", "right", "bottom", "top", "back", "front"};
    int direction = -1;
    for (int k = 0; k < 6; ++k)
      if (!strcmp (side, sides[k]))
        direction = k;
    if (direction < 0 || direction >= 2 * dimension) {
      basilisk_xml_error (
        error, line, "Unknown side or side outside the compiled dimension.");
      {
        basilisk_xml_case_destroy (&temporary);
        return 1;
      }
    }
    elff_xml_status status =
      elff_xml_node_unique_child (node, "Value", &value_node, error);
    if (!strcmp (type, "periodic")) {
      if (strcmp (field, "*") || status != ELFF_XML_NOT_FOUND) {
        basilisk_xml_error (
          error,
          line,
          "Periodicity requires field='*' and no Value; it affects the entire axis.");
        {
          basilisk_xml_case_destroy (&temporary);
          return 1;
        }
      }
      temporary.periodic_axis[direction / 2] = true;
      elff_xml_error_clear (error);
    } else {
      double value = 0.;
      if (status != ELFF_XML_OK ||
          elff_xml_node_text_double (value_node, &value, error) != ELFF_XML_OK) {
        basilisk_xml_case_destroy (&temporary);
        return 1;
      }
      BasiliskBCType kind;
      if (!strcmp (type, "dirichlet"))
        kind = BASILISK_BC_DIRICHLET;
      else if (!strcmp (type, "neumann"))
        kind = BASILISK_BC_NEUMANN;
      else {
        basilisk_xml_error (error, line, "Unknown boundary condition type.");
        {
          basilisk_xml_case_destroy (&temporary);
          return 1;
        }
      }
      if (basilisk_xml_add_bc (&temporary,
                               field,
                               (BasiliskBCDirection) direction,
                               kind,
                               value,
                               line,
                               error)) {
        basilisk_xml_case_destroy (&temporary);
        return 1;
      }
    }
  }
  // Existence checks are deferred until fields have been allocated.
  for (size_t j = 0; j < temporary.bc_count; ++j)
    if (temporary.periodic_axis[temporary.bcs[j].direction / 2]) {
      basilisk_xml_error (error,
                          temporary.bcs[j].source_line,
                          "Field boundary conflicts with a periodic axis.");
      {
        basilisk_xml_case_destroy (&temporary);
        return 1;
      }
    }
  free (c->bcs);
  c->bcs = temporary.bcs;
  c->bc_count = temporary.bc_count;
  memcpy (c->periodic_axis, temporary.periodic_axis, sizeof (c->periodic_axis));
  return 0;
}

static int generic_ibm_case_run (const BasiliskXMLCase* configuration) {
  elff_xml_error error = {0};
  if (!generic_case_agree (configuration != NULL))
    return 1;
  generic_case = *configuration;
  BasiliskXMLCase* c = &generic_case;
  bool valid =
    isfinite (c->time_max_dt) && c->time_max_dt > 0 && isfinite (c->time_end) &&
    c->time_end > 0 && isfinite (c->fluid_density) && c->fluid_density > 0 &&
    isfinite (c->dynamic_viscosity) && c->dynamic_viscosity >= 0 &&
    c->grid_min_level >= 0 && c->grid_max_level >= c->grid_min_level &&
    c->grid_max_level < 30 && c->grid_initial_level >= c->grid_min_level &&
    c->grid_initial_level <= c->grid_max_level && c->ib_level >= c->grid_min_level &&
    c->ib_level < 30 && c->force_iterations > 0 && isfinite (c->force_relaxation) &&
    c->force_relaxation > 0 && c->force_relaxation <= 1 &&
    generic_schedule_valid (c->fluid_output) && generic_schedule_valid (c->solid_output);
  for (int d = 0; d < dimension; ++d)
    valid = valid && isfinite (c->grid_size[d]) && c->grid_size[d] > 0 &&
            c->grid_size[d] == c->grid_size[0] && isfinite (c->grid_origin[d]) &&
            isfinite (c->force_density[d]);
  if (c->adapt_enable) {
    if (!c->adapt_fields || !c->adapt_tolerances ||
        c->adapt_tolerance_count != (size_t) list_len (c->adapt_fields))
      valid = false;
    else
      for (size_t k = 0; k < c->adapt_tolerance_count; ++k)
        if (!isfinite (c->adapt_tolerances[k]) || c->adapt_tolerances[k] <= 0)
          valid = false;
  }
#if !TREE
  if (c->adapt_enable)
    valid = false;
#endif
  if (valid)
    valid = !basilisk_xml_validate_boundaries (c, &error);
  if (!generic_case_agree (valid)) {
    if (!valid)
      fprintf (stderr,
               "Generic case configuration error: %s\n",
               error.message[0]
                 ? error.message
                 : "invalid grid, material, schedule or coupling settings");
    return 1;
  }
  size (c->grid_size[0]);
  origin (c->grid_origin[0], c->grid_origin[1], c->grid_origin[2]);
  if (c->periodic_axis[0])
    periodic (right);
#if dimension > 1
  if (c->periodic_axis[1])
    periodic (top);
#endif
#if dimension > 2
  if (c->periodic_axis[2])
    periodic (front);
#endif
  N = 1 << c->grid_initial_level;
  DT = c->time_max_dt;
  generic_end_time = c->time_end;
  // Disabled output must not impose extra timesteps.
  generic_fluid_start = c->fluid_output.enabled ? c->fluid_output.start : HUGE;
  generic_fluid_interval = c->fluid_output.enabled ? c->fluid_output.interval : 1.;
  generic_fluid_end = c->fluid_output.enabled ? c->fluid_output.end : HUGE;
  generic_solid_start = c->solid_output.enabled ? c->solid_output.start : HUGE;
  generic_solid_interval = c->solid_output.enabled ? c->solid_output.interval : 1.;
  generic_solid_end = c->solid_output.enabled ? c->solid_output.end : HUGE;
  mu = generic_mu;
  alpha = generic_alpha;
  rho = generic_rho;
  if (c->kernel)
    default_kernel = c->kernel;
  checkpoint_configuration.checkpoint_on_sim_time = c->checkpoint_interval > 0;
  checkpoint_configuration.checkpoint_on_sim_time_seconds = c->checkpoint_interval;
  // ib_force_relaxation = c->force_relaxation;
  // ib_richardson_iters = c->force_iterations;
  generic_case_status = 0;
  run ();
  return generic_case_status;
}

/** @brief Read an optional child without hiding malformed/duplicate-node errors. */
static int basilisk_xml_optional (elff_xml_node_t parent,
                                  const char* name,
                                  elff_xml_node_t* node,
                                  elff_xml_error* error) {
  elff_xml_status status = elff_xml_node_unique_child (parent, name, node, error);
  if (status == ELFF_XML_NOT_FOUND) {
    *node = NULL;
    elff_xml_error_clear (error);
    return 0;
  }
  return status != ELFF_XML_OK;
}

static int
basilisk_xml_enabled (elff_xml_node_t node, bool* enabled, elff_xml_error* error) {
  *enabled = node != NULL;
  if (!node)
    return 0;
  char text[16];
  elff_xml_status status =
    elff_xml_node_attribute_string (node, "enable", text, sizeof (text), NULL, error);
  if (status == ELFF_XML_NOT_FOUND) {
    elff_xml_error_clear (error);
    return 0;
  }
  if (status != ELFF_XML_OK)
    return 1;
  if (!strcmp (text, "true"))
    *enabled = true;
  else if (!strcmp (text, "false"))
    *enabled = false;
  else
    return basilisk_xml_error (error, 0, "enable must be true or false.");
  return 0;
}

static int basilisk_xml_schedule (elff_xml_node_t node,
                                  double end,
                                  BasiliskOutputSchedule* schedule,
                                  elff_xml_error* error) {
  schedule->start = 0.;
  schedule->interval = 1.;
  schedule->end = end;
  if (basilisk_xml_enabled (node, &schedule->enabled, error))
    return 1;
  if (!node)
    return 0;
  elff_xml_node_t time = NULL;
  if (basilisk_xml_optional (node, "SimulationTime", &time, error))
    return 1;
  if (time &&
      (elff_xml_node_attribute_double (time, "step", &schedule->interval, error) !=
         ELFF_XML_OK ||
       elff_xml_node_attribute_double_or (time, "end", end, &schedule->end, error) !=
         ELFF_XML_OK))
    return 1;
  // Retain the historical first output at one interval, not at t = 0.
  schedule->start = schedule->interval;
  if (!generic_schedule_valid (*schedule))
    return basilisk_xml_error (error, 0, "Invalid simulation-time output schedule.");
  return 0;
}

/* Read names now; binding is delayed until solver defaults allocate IB fields. */
static int basilisk_xml_output_fields (elff_xml_node_t node,
                                       bool solid,
                                       BasiliskXMLField** fields,
                                       size_t* count,
                                       elff_xml_error* error) {
  if (!node)
    return 0;
  size_t n = 0;
  if (elff_xml_node_element_child_count (node, &n, error) != ELFF_XML_OK)
    return 1;
  if (n > (size_t) -1 / sizeof (BasiliskXMLField))
    return basilisk_xml_error (error, 0, "Too many output fields.");
  *fields = calloc (n ? n : 1, sizeof (BasiliskXMLField));
  if (!*fields)
    return basilisk_xml_error (error, 0, "Cannot allocate output fields.");
  for (size_t j = 0; j < n; ++j) {
    elff_xml_node_t child = NULL;
    const char* tag = NULL;
    if (elff_xml_node_element_child_at (node, j, &child, error) != ELFF_XML_OK ||
        elff_xml_node_name (child, &tag, error) != ELFF_XML_OK)
      return 1;
    if (!strcmp (tag, "SimulationTime"))
      continue;
    bool vector_field = !strcmp (tag, solid ? "IBVectorField" : "VectorField");
    if (!vector_field && strcmp (tag, solid ? "IBScalarField" : "ScalarField"))
      return basilisk_xml_error (error, 0, "Unsupported output field element.");
    BasiliskXMLField* field = &(*fields)[(*count)++];
    field->vector_field = vector_field;
    elff_xml_node_line (child, &field->source_line, NULL);
    if (elff_xml_node_attribute_string (
          child, "name", field->name, sizeof (field->name), NULL, error) != ELFF_XML_OK)
      return 1;
    for (size_t k = 0; k + 1 < *count; ++k)
      if (!strcmp ((*fields)[k].name, field->name))
        return basilisk_xml_error (
          error, field->source_line, "Duplicate output field name.");
  }
  return 0;
}

static int basilisk_xml_bind_output_fields (BasiliskXMLCase* c, elff_xml_error* error) {
  if (!c->owns_field_lists)
    return 0;
  for (size_t j = 0; j < c->fluid_field_count; ++j) {
    BasiliskXMLField* field = &c->fluid_fields[j];
    if (field->vector_field) {
      vector v = lookup_vector (field->name);
      if (v.x.i < 0)
        return basilisk_xml_error (
          error, field->source_line, "Unknown fluid output vector.");
      c->fluid_vectors = vectors_append (c->fluid_vectors, v);
    } else {
      scalar s = lookup_field (field->name);
      if (s.i < 0)
        return basilisk_xml_error (
          error, field->source_line, "Unknown fluid output scalar.");
      c->fluid_scalars = list_append (c->fluid_scalars, s);
    }
  }
  for (size_t j = 0; j < c->solid_field_count; ++j) {
    BasiliskXMLField* field = &c->solid_fields[j];
    char name[132];
    snprintf (name, sizeof (name), "%s%s", field->name, field->vector_field ? ".x" : "");
    IBscalar s = {-1};
    for (size_t k = 0; k < _ibattribute_len; ++k)
      if (_ibattribute[k].name && !strcmp (_ibattribute[k].name, name)) {
        s.i = k;
        break;
      }
    if (s.i < 0)
      return basilisk_xml_error (error, field->source_line, "Unknown IB output field.");
    if (field->vector_field) {
      IBvector v = _ibattribute[s.i].v;
      if (v.x.i != s.i)
        return basilisk_xml_error (
          error, field->source_line, "IB field is not a vector.");
      c->solid_vectors = ibvectors_append (c->solid_vectors, v);
    } else
      c->solid_scalars = iblist_append (c->solid_scalars, s);
  }
  return 0;
}

static int basilisk_xml_case_read (elff_xml_document_t document,
                                   const char* expected_case,
                                   BasiliskXMLCase* configuration,
                                   elff_xml_error* error) {
  BasiliskXMLCase parsed;
  basilisk_xml_case_defaults (&parsed);
  BasiliskXMLCase* c = &parsed;
  c->owns_field_lists = true;
  elff_xml_node_t root = NULL, b = NULL, time = NULL, fluid = NULL, material = NULL;
  elff_xml_node_t grid = NULL, node = NULL, adaptation = NULL, ib = NULL, output = NULL;
  int version = 0;
#define XML_READ(call)                                                                   \
  do {                                                                                   \
    if ((call) != ELFF_XML_OK) {                                                         \
      basilisk_xml_case_destroy (&parsed);                                               \
      return 1;                                                                          \
    }                                                                                    \
  } while (0)
#define XML_CHECK(call)                                                                  \
  do {                                                                                   \
    if (call) {                                                                          \
      basilisk_xml_case_destroy (&parsed);                                               \
      return 1;                                                                          \
    }                                                                                    \
  } while (0)
  XML_READ (elff_xml_case_validate (document, error));
  XML_READ (elff_xml_document_root (document, "Case", &root, error));
  XML_READ (elff_xml_node_attribute_int (root, "formatVersion", &version, error));
  if (version != 2) {
    basilisk_xml_error (error, 0, "This executable requires XML formatVersion=2.");
    {
      basilisk_xml_case_destroy (&parsed);
      return 1;
    }
  }
  XML_CHECK (basilisk_xml_validate_backend (document, error));
  XML_READ (basilisk_xml_read_named_values (document, c, error));
  XML_READ (elff_xml_node_unique_child (root, "Basilisk", &b, error));
  XML_READ (elff_xml_node_attribute_string (
    b, "caseType", c->case_name, sizeof (c->case_name), NULL, error));
  if (!expected_case || strcmp (c->case_name, expected_case)) {
    size_t line = 0;
    elff_xml_node_line (b, &line, NULL);
    basilisk_xml_error (error, line, "XML caseType does not match this executable.");
    {
      basilisk_xml_case_destroy (&parsed);
      return 1;
    }
  }
  XML_READ (elff_xml_node_unique_child (b, "Time", &time, error));
  XML_READ (elff_xml_node_attribute_double (time, "step", &c->time_max_dt, error));
  XML_READ (elff_xml_node_attribute_double (time, "end", &c->time_end, error));
  XML_READ (elff_xml_node_unique_child (b, "Fluid", &fluid, error));
  XML_READ (elff_xml_node_unique_child (fluid, "FluidMaterial", &material, error));
  XML_READ (
    elff_xml_node_attribute_double (material, "density", &c->fluid_density, error));
  XML_READ (elff_xml_node_attribute_double (
    material, "dynamicViscosity", &c->dynamic_viscosity, error));
  XML_CHECK (basilisk_xml_read_boundaries (fluid, c, error));
  XML_READ (elff_xml_node_unique_child (b, "Grid", &grid, error));
  XML_READ (elff_xml_node_unique_child (grid, "Origin", &node, error));
  XML_READ (elff_xml_node_attribute_xyz (node, 0., c->grid_origin, error));
  XML_READ (elff_xml_node_unique_child (grid, "Length", &node, error));
  XML_READ (elff_xml_node_attribute_xyz (node, 1., c->grid_size, error));
#if dimension == 3
  XML_READ (elff_xml_node_attribute_double (node, "z", &c->grid_size[2], error));
#endif
  XML_READ (elff_xml_node_unique_child (grid, "Depth", &node, error));
  double depth = 0.;
  XML_READ (elff_xml_node_text_double (node, &depth, error));
  if (depth < 0 || depth >= 30 || floor (depth) != depth) {
    basilisk_xml_error (error, 0, "Depth must be an integer in [0,30).");
    {
      basilisk_xml_case_destroy (&parsed);
      return 1;
    }
  }
  c->grid_initial_level = c->grid_min_level = c->grid_max_level = (int) depth;
  XML_CHECK (basilisk_xml_optional (grid, "AdaptWavelet", &adaptation, error));
  if (adaptation) {
    c->adapt_enable = true;
    XML_READ (
      elff_xml_node_attribute_int (adaptation, "minLevel", &c->grid_min_level, error));
    XML_READ (
      elff_xml_node_attribute_int (adaptation, "maxLevel", &c->grid_max_level, error));
    size_t count = 0;
    XML_READ (elff_xml_node_element_child_count (adaptation, &count, error));
    if (count > (size_t) -1 / (dimension * sizeof (double))) {
      basilisk_xml_error (error, 0, "Too many adaptation fields.");
      {
        basilisk_xml_case_destroy (&parsed);
        return 1;
      }
    }
    c->adapt_tolerances = calloc (count ? count * dimension : 1, sizeof (double));
    if (!c->adapt_tolerances) {
      basilisk_xml_error (error, 0, "Cannot allocate tolerances.");
      {
        basilisk_xml_case_destroy (&parsed);
        return 1;
      }
    }
    for (size_t j = 0; j < count; ++j) {
      const char* tag = NULL;
      char name[128];
      double tolerance = 0.;
      XML_READ (elff_xml_node_element_child_at (adaptation, j, &node, error));
      XML_READ (elff_xml_node_name (node, &tag, error));
      XML_READ (
        elff_xml_node_attribute_string (node, "name", name, sizeof (name), NULL, error));
      XML_READ (elff_xml_node_attribute_double (node, "tolerance", &tolerance, error));
      if (!strcmp (tag, "VectorField")) {
        vector v = lookup_vector (name);
        if (v.x.i < 0) {
          basilisk_xml_error (error, 0, "Unknown adaptation vector.");
          {
            basilisk_xml_case_destroy (&parsed);
            return 1;
          }
        }
        foreach_dimension () {
          c->adapt_fields = list_append (c->adapt_fields, v.x);
          c->adapt_tolerances[c->adapt_tolerance_count++] = tolerance;
        }
      } else {
        scalar s = lookup_field (name);
        if (s.i < 0) {
          basilisk_xml_error (error, 0, "Unknown adaptation scalar.");
          {
            basilisk_xml_case_destroy (&parsed);
            return 1;
          }
        }
        c->adapt_fields = list_append (c->adapt_fields, s);
        c->adapt_tolerances[c->adapt_tolerance_count++] = tolerance;
      }
    }
  }
  c->ib_level = c->grid_max_level;
  XML_CHECK (basilisk_xml_optional (b, "ImmersedBoundary", &ib, error));
  if (ib) {
    elff_xml_status status =
      elff_xml_node_attribute_int (ib, "level", &c->ib_level, error);
    if (status == ELFF_XML_NOT_FOUND)
      elff_xml_error_clear (error);
    else if (status != ELFF_XML_OK) {
      basilisk_xml_case_destroy (&parsed);
      return 1;
    }
    XML_CHECK (basilisk_xml_optional (ib, "UniformForceDensity", &node, error));
    if (node)
      XML_READ (elff_xml_node_attribute_xyz (node, 0., c->force_density, error));
    XML_CHECK (basilisk_xml_optional (ib, "ForceCoupled", &node, error));
    if (node) {
      elff_xml_node_t control = NULL;
      XML_CHECK (basilisk_xml_optional (node, "Iterations", &control, error));
      if (control) {
        double iterations = 0.;
        XML_READ (elff_xml_node_text_double (control, &iterations, error));
        if (iterations < 1 || iterations > INT_MAX || floor (iterations) != iterations) {
          basilisk_xml_error (error, 0, "Invalid force iteration count.");
          {
            basilisk_xml_case_destroy (&parsed);
            return 1;
          }
        }
        c->force_iterations = (int) iterations;
      }
      XML_CHECK (basilisk_xml_optional (node, "Relaxation", &control, error));
      if (control)
        XML_READ (elff_xml_node_text_double (control, &c->force_relaxation, error));
    }
    XML_CHECK (basilisk_xml_optional (ib, "ImmersedBoundaryKernel", &node, error));
    if (node) {
      char name[64];
      XML_READ (
        elff_xml_node_attribute_string (node, "type", name, sizeof (name), NULL, error));
      if (!strcmp (name, "Roma"))
        c->kernel = roma_three_point_kernel_1d;
      else if (!strcmp (name, "Cosine"))
        c->kernel = peskin_cosine_kernel_1d;
      else {
        basilisk_xml_error (error, 0, "Kernel must be Roma or Cosine.");
        {
          basilisk_xml_case_destroy (&parsed);
          return 1;
        }
      }
    }
  }
  XML_CHECK (basilisk_xml_optional (b, "Restart", &node, error));
  XML_CHECK (basilisk_xml_enabled (node, &c->restore_enable, error));
  XML_READ (elff_xml_node_unique_child (b, "Output", &output, error));
  XML_READ (elff_xml_node_attribute_string (
    output, "directory", c->output_directory, sizeof (c->output_directory), NULL, error));
  XML_CHECK (basilisk_xml_optional (output, "FluidOutput", &node, error));
  XML_CHECK (basilisk_xml_schedule (node, c->time_end, &c->fluid_output, error));
  XML_CHECK (basilisk_xml_output_fields (
    node, false, &c->fluid_fields, &c->fluid_field_count, error));
  XML_CHECK (basilisk_xml_optional (output, "SolidOutput", &node, error));
  XML_CHECK (basilisk_xml_schedule (node, c->time_end, &c->solid_output, error));
  XML_CHECK (basilisk_xml_output_fields (
    node, true, &c->solid_fields, &c->solid_field_count, error));
  XML_CHECK (basilisk_xml_optional (output, "Checkpoint", &node, error));
  XML_CHECK (basilisk_xml_enabled (node, &c->checkpoint_enable, error));
  if (node) {
    elff_xml_node_t schedule = NULL;
    XML_CHECK (basilisk_xml_optional (node, "SimulationTime", &schedule, error));
    if (schedule) {
      XML_READ (elff_xml_node_attribute_double (
        schedule, "step", &c->checkpoint_interval, error));
      XML_READ (elff_xml_node_attribute_double_or (
        schedule, "end", c->time_end, &c->checkpoint_end, error));
      if (c->checkpoint_interval <= 0 || c->checkpoint_end < 0) {
        basilisk_xml_error (error, 0, "Invalid checkpoint schedule.");
        {
          basilisk_xml_case_destroy (&parsed);
          return 1;
        }
      }
    }
    size_t count = 0;
    XML_READ (elff_xml_node_element_child_count (node, &count, error));
    if (count > (schedule ? 1u : 0u)) {
      basilisk_xml_error (
        error, 0, "Additional checkpoint field selection is not implemented.");
      {
        basilisk_xml_case_destroy (&parsed);
        return 1;
      }
    }
  }
  basilisk_xml_case_destroy (configuration);
  *configuration = parsed;
  return 0;
#undef XML_READ
#undef XML_CHECK
}
