#include "library/ibm/IBMeshManager.h"
#include "library/ibm/navier-stokes/centered-mdf.h"
#include "library/elff/elff.h"
#include "library/ibm/IBOutput.h"
#include "library/io/output-vtk.h"
#include "elff/io/C/xmlCase.h"
#include "elff/io/DTD/embeddedFlappingFibreDtd.h"

#include <math.h>
#include <string.h>

/* 
  The default simulation parameters correspond to Figure 13 of https://doi.org/10.1016/j.jcp.2007.07.002 
*/

typedef struct
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
  int owner_rank;
  int has_uniform_force_density;
  double uniform_force_density[3];
  char output_directory[256];
  double output_interval;
} FlappingFibreConfiguration;

static elff_xml_status
flapping_fibre_error (elff_xml_error * error, const char * message)
{
  elff_xml_error_set (error, 0, 0, message);
  return ELFF_XML_INVALID_VALUE;
}

static elff_xml_status
flapping_fibre_close_document (elff_xml_document_t document,
                               elff_xml_status status)
{
  elff_xml_document_destroy (document);
  return status;
}

static elff_xml_status
flapping_fibre_case_load (const char * path,
                          FlappingFibreConfiguration * configuration,
                          elff_xml_error * error)
{
  elff_xml_document_t document = NULL;
  elff_xml_node_t basilisk = NULL, domain = NULL, time = NULL;
  elff_xml_node_t fluid = NULL, inlet = NULL, grid = NULL, immersed = NULL;
  elff_xml_node_t output = NULL, force = NULL;
  elff_xml_status status;

  if (configuration != NULL)
    memset (configuration, 0, sizeof (*configuration));
  if (path == NULL || configuration == NULL)
    return flapping_fibre_error (error,
                                 "Invalid flapping-fibre case argument.");

  status = elff_xml_case_load_section_file (
    path,
    "Case",
    "Basilisk",
    (const char *) elff_embedded_flapping_fibre_v1_dtd,
    elff_embedded_flapping_fibre_v1_dtd_size,
    &document,
    &basilisk,
    error);
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
         basilisk, "ImmersedBoundary", &immersed, error)) != ELFF_XML_OK ||
      (status = elff_xml_node_unique_child (
         basilisk, "Output", &output, error)) != ELFF_XML_OK)
    return flapping_fibre_close_document (document, status);

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
      (status = elff_xml_node_attribute_int (
         immersed,
         "level",
         &configuration->immersed_boundary_level,
         error)) != ELFF_XML_OK ||
      (status = elff_xml_node_attribute_int (
         immersed, "ownerRank", &configuration->owner_rank, error)) != ELFF_XML_OK ||
      (status = elff_xml_node_attribute_double (
         output, "interval", &configuration->output_interval, error)) != ELFF_XML_OK)
    return flapping_fibre_close_document (document, status);

  status = elff_xml_node_attribute_string (
    output,
    "directory",
    configuration->output_directory,
    sizeof (configuration->output_directory),
    NULL,
    error);
  if (status != ELFF_XML_OK)
    return flapping_fibre_close_document (document, status);

  status = elff_xml_node_unique_child (
    immersed, "UniformForceDensity", &force, error);
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
      return flapping_fibre_close_document (document, status);
  }
  else
    return flapping_fibre_close_document (document, status);

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
      configuration->owner_rank != 0 ||
      !(configuration->output_interval > 0.))
    return flapping_fibre_close_document (
      document,
      flapping_fibre_error (
        error, "Invalid numeric value in <Basilisk> configuration."));

  return flapping_fibre_close_document (document, ELFF_XML_OK);
}

static FlappingFibreConfiguration case_configuration = {0};
static char* base_path = NULL;
static double U0 = 0.;
static double fluid_density = 1.;
static double fluid_viscosity = 0.;
static double output_interval = 0.05;
static double end_time = 0.;
static int maxlevel = 0;
static int minlevel = 0;
static int ibmlevel = 0;

/* Additional fields */
face vector muv[], alphav[];
scalar rhov[];

/* Boundary conditions */

u.n[left] = dirichlet(U0);
u.t[left] = dirichlet(0.);
p[left] = neumann(0.);
pf[left] = neumann(0.);

u.n[right] = neumann(0.);
p[right] = dirichlet(0.);
pf[right] = dirichlet(0.);


int
main(int argc, char ** argv)
{
  const char* case_path = argc > 1 ? argv[1] : "flapping_fibre.xml";
  elff_xml_error basilisk_error = {0};
  if (flapping_fibre_case_load (
        case_path, &case_configuration, &basilisk_error) != ELFF_XML_OK) {
    fprintf (stderr, "Basilisk case error: %s\n", basilisk_error.message);
    return 1;
  }
  ib_runtime_error runtime_error = {0};
  if (elff_init (case_path, &runtime_error) != 0) {
    fprintf (stderr,
             "ELFF model %zu (line %ld): %s\n",
             runtime_error.model_index,
             runtime_error.line,
             runtime_error.message);
    return 1;
  }

  /* Setting relevant parameters for basilisk */
  L0 = case_configuration.domain_size;
  origin(case_configuration.origin[0], case_configuration.origin[1]);
  U0 = case_configuration.inlet_velocity[0];
  fluid_density = case_configuration.fluid_density;
  fluid_viscosity = case_configuration.dynamic_viscosity;
  minlevel = case_configuration.min_level;
  maxlevel = case_configuration.max_level;
  ibmlevel = case_configuration.immersed_boundary_level;
  base_path = case_configuration.output_directory;
  output_interval = case_configuration.output_interval;
  end_time = case_configuration.end_time;
#if TREE
  N = 1 << minlevel;
#else 
  N = 1 << ibmlevel;
#endif   
  DT = case_configuration.time_step;
  mu = muv;
  alpha = alphav;
  rho = rhov;

  run();
}

/* Reynolds control */
event
properties(i++)
{
  foreach_face() {
    muv.x[] = fm.x[] * fluid_viscosity;
    alphav.x[] = fm.x[] / fluid_density;
  }
  foreach()
    rhov[] = cm[] * fluid_density;
}

/* Beam setup */
event
init(i = 0)
{
  ib_model_view_array_t models = elff_models ();
  assert (models.count > 0);
  for (size_t model_index = 0; model_index < models.count; ++model_index) {
    int mesh_index = ibmeshmanager_add_mesh ();
    assert (mesh_index == (int) model_index);
    IBMeshModel adapter = elff_model_from_view (models.data[model_index]);
    assert (adapter.ctx);
    ibmeshmanager_set_model (mesh_index, adapter);
  }

  foreach_ibnode_per_ibmesh()
  {
    mesh->depth = ibmlevel;
    node->depth = ibmlevel;
    if (case_configuration.has_uniform_force_density) {
      ibval(gravity.x) = case_configuration.uniform_force_density[0];
      ibval(gravity.y) = case_configuration.uniform_force_density[1];
#if dimension > 2
      ibval(gravity.z) = case_configuration.uniform_force_density[2];
#endif
    }
  }

  if (!restore_handler(base_path)) {
    foreach ()
      u.x[] = U0;
#if TREE
    adapt_wavelet_ibm(NULL, NULL, 0, 1, all, true);
#endif
  } else {
  }
}

event
logfile(i++) {
  if (pid() == 0)
    fprintf(stderr, "%d %g\n", i, t);
}

event
statsfile(i++)
{
  double x_tip = 0.;
  double y_tip = 0.;
  {
    IBNode* node = ibmm.pool.active.ptrs[0];
    x_tip = ibval(npos.x);
    y_tip = ibval(npos.y);
  }

  if (pid() == 0) {
    FILE *fp = NULL;
    if (!base_path) {
      fprintf(stderr, "warning: base_path is NULL; skipping tip output\n");
      return 0;
    }
    create_path(base_path);
    char fname[4096];
    snprintf(fname, sizeof(fname), "%s/tip.csv", base_path);
    fp = fopen(fname, "a");
    if (!fp) {
      fprintf(stderr, "warning: failed to open %s for append\n", fname);
      return 0;
    }
    if (i==0) {
      fprintf(fp, "i,t,x_tip,y_tip\n");
    }
    fprintf(fp, "%d,%g,%g,%g\n", i, t, x_tip, y_tip);
    fclose(fp);
  }
}


event
output(t += output_interval; t <= end_time)
{
  scalar omega[];
  vorticity(u, omega);
#if TREE
  output_hdf_htg({omega,p}, {u,ibmf}, base_path);
#else
  output_hdf_imagedata({omega,p}, {u,ibmf}, base_path);
#endif
  output_hdf_pd(NULL, (IBvector[]){eulvel, nforce, nvel, {{-1}}}, base_path);
}

#if TREE
event
adapt(i++)
{
  adapt_wavelet_ibm({ u }, (double[]){ 5e-3, 5e-3 }, maxlevel, minlevel);
}
#endif

event
checkpoint_event(i++, last)
{
  return checkpoint_handler(t, i, base_path);
}

event
elff_cleanup(i = end, last)
{
  ibmeshmanager_free ();
  elff_deinit ();
}
