
#include "library/ibm/navier-stokes/centered-mdf.h"
#include "library/io/xml/xml.h"

// #include "flapping_fibre_parser.h"

typedef struct {
  BasiliskXMLCase basilisk;
  double initial_velocity;
} FlappingFibreCase;

static int flapping_fibre_case_init(const char *path,
                                    FlappingFibreCase *configuration) {
  basilisk_xml_case_defaults(&configuration->basilisk);
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
        basilisk_xml_coefficient(&configuration->basilisk, "inletSpeed",
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
