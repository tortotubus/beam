
#include "library/ibm/navier-stokes/centered-mdf.h"
#include "library/io/xml/xml.h"

// #include "flapping_fibre_parser.h"

typedef struct {
  BasiliskXMLCase basilisk;
  double initial_velocity;
} FlappingFibreCase;

/** @brief Default nondimensional flapping domain, inlet and fluid controls. */
static int flapping_fibre_defaults (BasiliskXMLCase* c) {
  basilisk_xml_case_defaults(c);
  for (int d = 0; d < dimension; ++d) { c->grid_size[d] = 16.; c->grid_origin[d] = -8.; }
  c->grid_origin[0] = -2.;
  c->grid_initial_level = c->grid_min_level = 6;
  c->grid_max_level = 10; c->ib_level = 11;
  c->time_max_dt = 0.0003; c->time_end = 50.;
  c->fluid_density = 1.; c->dynamic_viscosity = 0.005;
  c->force_density[0] = 0.75; c->force_relaxation = 1.;
  c->fluid_output.interval = c->solid_output.interval = 0.05;
  strcpy(c->output_directory, "flapping_fibre_output");
  const BasiliskXMLField fields[] = {{"omega", false, 0}, {"p", false, 0},
    {"u", true, 0}, {"ibmf", true, 0}};
  if (basilisk_xml_default_fields(c, fields, 4, 0.005)) return 1;
  elff_xml_error error = {0};
  if (basilisk_xml_add_bc(c, "u.x", BASILISK_BC_LEFT, BASILISK_BC_DIRICHLET, 1., 0, &error) ||
      basilisk_xml_add_bc(c, "u.y", BASILISK_BC_LEFT, BASILISK_BC_DIRICHLET, 0., 0, &error) ||
      basilisk_xml_add_bc(c, "p", BASILISK_BC_LEFT, BASILISK_BC_NEUMANN, 0., 0, &error) ||
      basilisk_xml_add_bc(c, "pf", BASILISK_BC_LEFT, BASILISK_BC_NEUMANN, 0., 0, &error) ||
      basilisk_xml_add_bc(c, "u.x", BASILISK_BC_RIGHT, BASILISK_BC_NEUMANN, 0., 0, &error) ||
      basilisk_xml_add_bc(c, "p", BASILISK_BC_RIGHT, BASILISK_BC_DIRICHLET, 0., 0, &error) ||
      basilisk_xml_add_bc(c, "pf", BASILISK_BC_RIGHT, BASILISK_BC_DIRICHLET, 0., 0, &error)) {
    basilisk_xml_case_destroy(c); return 1;
  }
  return 0;
}

static int flapping_fibre_case_init(const char *path,
                                    FlappingFibreCase *configuration) {
  int defaults_ok = !flapping_fibre_defaults(&configuration->basilisk);
  if (!generic_case_agree(defaults_ok)) {
    basilisk_xml_case_destroy(&configuration->basilisk); return 1;
  }
  configuration->initial_velocity = 1.;

  elff_xml_document_t document = NULL;
  elff_xml_error error = {0};

  int success = elff_xml_case_load_rendered_file(path, pid() == 0, &document,
                                                 &error) == ELFF_XML_OK;
  if (success)
    success = !basilisk_xml_case_read(document, "flapping_fibre",
                                      &configuration->basilisk, &error);
  if (success) {
    elff_xml_status status =
        basilisk_xml_parameter(&configuration->basilisk, "inletSpeed",
                                 &configuration->initial_velocity, &error);
    success = status == ELFF_XML_OK || status == ELFF_XML_NOT_FOUND;
    if (status == ELFF_XML_NOT_FOUND)
      elff_xml_error_clear(&error);
  }

  if (!generic_case_agree(success)) {
    if (!success)
      fprintf(stderr, "flapping_fibre XML error (line %zu): %s\n", error.line,
              error.message);
    elff_xml_document_destroy(document);
    basilisk_xml_case_destroy(&configuration->basilisk);
    return 1;
  }

  model_runtime_error runtime_error = {0};
  int status = elff_init_document(document, &runtime_error);
  elff_xml_document_destroy(document);

  if (status) {
    fprintf(stderr, "ELFF model %zu (line %ld): %s\n",
            runtime_error.model_index, runtime_error.line,
            runtime_error.message);
    basilisk_xml_case_destroy(&configuration->basilisk);
    return status;
  }

  ib_model_view_array_t models = elff_models();

  if (!models.count || models.count > INT_MAX) {
    fprintf(stderr, "flapping_fibre requires a nonempty model population "
                    "fitting mesh indices.\n");
    elff_deinit();
    basilisk_xml_case_destroy(&configuration->basilisk);
    return 1;
  }

  return 0;
}

static FlappingFibreCase case_configuration = {0};
static char *base_path = NULL;
scalar omega[];

static void flapping_initialize_fluid(void) {
  foreach ()
    u.x[] = case_configuration.initial_velocity;
}

static void flapping_prepare_output(void) { vorticity(u, omega); }

int main(int argc, char **argv) {
  const char *path = argc > 1 ? argv[1] : "flapping_fibre.xml";
  if (flapping_fibre_case_init(path, &case_configuration))
    return 1;

  BasiliskXMLCase *c = &case_configuration.basilisk;
  base_path = c->output_directory;

  c->initialize_fluid = flapping_initialize_fluid;
  c->prepare_fluid_output = flapping_prepare_output;
  ib_force_relaxation = c->force_relaxation;
  ib_richardson_iters = c->force_iterations;  

  int status = generic_ibm_case_run(c);


  elff_deinit();
  basilisk_xml_case_destroy(c);
  return status;
}

event logfile(i++) {
  if (pid() == 0)
    fprintf(stderr, "%d %g\n", i, t);
}

event statsfile(i++) {
  double x_tip = 0.;
  double y_tip = 0.;
  {
    foreach_ibnode() {
      if (node->mesh_gid == 0 && node->node_lid == 0 && node->pid == pid()) {
        x_tip = ibval(npos.x);
        y_tip = ibval(npos.y);
      }
    }
#if _MPI
    MPI_Allreduce(MPI_IN_PLACE, &x_tip, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &y_tip, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
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
    if (i == 0) {
      fprintf(fp, "i,t,x_tip,y_tip\n");
    }
    fprintf(fp, "%d,%g,%g,%g\n", i, t, x_tip, y_tip);
    fclose(fp);
  }
}
