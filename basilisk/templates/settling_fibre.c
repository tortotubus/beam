/** Shared XML-driven settling experiment for one or many fibres.
 * Banaei parameters and dimensional coefficients are evaluated in the case XML.
 */
#ifndef IBMESHMANAGER_POOL_SIZE_BYTES
#define IBMESHMANAGER_POOL_SIZE_BYTES (1 << 24)
#endif
#include "library/ibm/navier-stokes/centered-mdf.h"
#include "library/io/xml/xml.h"
#include "lambda2.h"
#if dimension != 3
#error "settling_fibre requires grid/octree.h"
#endif
#include "settling_fibre_parser.h"

static SettlingFibreCase settling_configuration;
static int stats_interval = 1;
static char *base_path;
scalar l2[], omega_z[];

/** @brief Refresh derived fields immediately before fluid output. */
static void settling_fibre_prepare_output () {
  lambda2 (u, l2);
  vorticity (u, omega_z);
}

static void settling_fibre_stats (int iteration, double time) {
  enum {
    Y_FIRST_SUM,
    Y_FIRST_COUNT,
    Y_LAST_SUM,
    Y_LAST_COUNT,
    FORCE_X,
    FORCE_Y,
    FORCE_Z,
    EXTERNAL_POWER,
    STATS_STRIDE
  };
  size_t count = (size_t) ibmm.nm * STATS_STRIDE;
  double *minimum_y = malloc ((size_t) ibmm.nm * sizeof (double));
  double *values = calloc (count, sizeof (double));
  if (!minimum_y || !values) {
    free (minimum_y);
    free (values);
    fprintf (stderr, "Unable to allocate fibre statistics.\n");
    exit (1);
  }
  for (int m = 0; m < ibmm.nm; ++m)
    minimum_y[m] = HUGE_VAL;
  // Count each MPI-owned node once. At this point in the event sequence,
  // nforce is the total external force density used to advance the fibre:
  // the previous IBM reaction plus the imposed body-force density. nweight
  // supplies the nodal measure for force and power integration.
  foreach_ibnode_per_ibmesh () {
#if _MPI
    if (node->pid != pid ()) continue;
#endif
    double *v = values + STATS_STRIDE*mesh_id;
    double y = ibval (npos.y);
    double w = ibval (nweight);
    if (y < minimum_y[mesh_id]) minimum_y[mesh_id] = y;
    if (node->node_lid == 0) {
      v[Y_FIRST_SUM] += y;
      v[Y_FIRST_COUNT] += 1.;
    }
    if (node->node_lid == mesh->nnode_global - 1) {
      v[Y_LAST_SUM] += y;
      v[Y_LAST_COUNT] += 1.;
    }
    v[FORCE_X] += w*ibval (nforce.x);
    v[FORCE_Y] += w*ibval (nforce.y);
    v[FORCE_Z] += w*ibval (nforce.z);
    v[EXTERNAL_POWER] += w*(ibval (nforce.x)*ibval (nvel.x) +
                            ibval (nforce.y)*ibval (nvel.y) +
                            ibval (nforce.z)*ibval (nvel.z));
  }
#if _MPI
  MPI_Allreduce (MPI_IN_PLACE, minimum_y, ibmm.nm, MPI_DOUBLE, MPI_MIN,
                 MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, values, (int)count, MPI_DOUBLE, MPI_SUM,
                 MPI_COMM_WORLD);
#endif
  if (pid () == 0) {
    create_path (base_path);
    char filename[4096];
    snprintf (filename, sizeof(filename), "%s/fibre-diagnostics.csv", base_path);
    FILE *fp = fopen (filename, iteration == 0 ? "w" : "a");
    if (!fp) { perror (filename); exit (1); }
    if (iteration == 0)
      fprintf (fp, "i,t,mesh_id,y_min,y_first_node,y_last_node,external_force_x,external_force_y,external_force_z,external_power\n");
    for (int m = 0; m < ibmm.nm; ++m) {
      double *v = values + STATS_STRIDE*m;
      double y_first = v[Y_FIRST_COUNT] > 0. ?
        v[Y_FIRST_SUM]/v[Y_FIRST_COUNT] : NAN;
      double y_last = v[Y_LAST_COUNT] > 0. ?
        v[Y_LAST_SUM]/v[Y_LAST_COUNT] : NAN;
      double y_min = minimum_y[m] < HUGE_VAL ? minimum_y[m] : NAN;
      fprintf (fp, "%d,%.17g,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
               iteration, time, m, y_min, y_first, y_last,
               v[FORCE_X], v[FORCE_Y], v[FORCE_Z], v[EXTERNAL_POWER]);
    }
    fclose (fp);
    fprintf (stderr, "settling i=%d t=%g fibres=%d\n", iteration, time, ibmm.nm);
  }
  free (minimum_y);
  free (values);
}

int main (int argc, char **argv) {
  if (argc != 2) {
    if (pid() == 0) fprintf (stderr, "Usage: settling_fibre case.xml\n");
    return 1;
  }
  if (settling_fibre_case_init (argv[1], &settling_configuration)) return 1;
  BasiliskXMLCase *configuration = &settling_configuration.basilisk;
  stats_interval = settling_configuration.stats_interval;
  base_path = configuration->output_directory;
  configuration->prepare_fluid_output = settling_fibre_prepare_output;
  ib_model_view_array_t models = elff_models ();
  if (pid() == 0)
    fprintf (stderr, "settling fibres=%zu periodic rho=%g mu=%g gravity_y=%g\n",
             models.count, configuration->fluid_density, configuration->dynamic_viscosity,
             configuration->force_density[1]);
  ib_force_relaxation = configuration->force_relaxation;
  ib_richardson_iters = configuration->force_iterations;
  int status = generic_ibm_case_run (configuration);
  elff_deinit ();
  basilisk_xml_case_destroy (configuration);
  return status;
}

event logfile (i += stats_interval; t <= generic_case.time_end) {
  settling_fibre_stats (i, t);
}
