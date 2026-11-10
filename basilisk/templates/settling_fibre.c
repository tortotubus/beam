/** Shared XML-driven settling experiment for one or many fibres.
 * Banaei parameters and dimensional coefficients are evaluated in the case XML.
 *
 * The MDF exchange diagnostic brackets interface_force_richardson, before
 * beta viscosity or projection can alter the result.  It reports
 *
 *   W_fluid = K_fluid(after MDF) - K_fluid(before MDF),
 *   W_fibre = dt sum_a w_a F_a^IBM . v_a,
 *   R_IBM   = W_fluid + W_fibre.
 *
 * R_IBM is the model-independent discrete interface-energy residual.  The
 * imposed UniformForceDensity is the only external source for the combined
 * fluid/fibre system; its power and step work are reported separately.  This
 * tests coupling work, not the fibre's kinetic plus elastic energy, which is
 * model-specific and is not exposed by every ELFF structural model.
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

static double fluid_kinetic_before_mdf = NAN;
static double fluid_ibm_work = NAN;
static double fibre_ibm_power = NAN;
static double external_power = NAN;
static double energy_step_dt = NAN;

/** @brief Cell-centred kinetic energy with the same rho and dv as the MDF update. */
static double settling_fluid_kinetic_energy () {
  double kinetic = 0.;
  foreach (reduction(+:kinetic)) {
    double speed2 = 0.;
    foreach_dimension () speed2 += sq (u.x[]);
    kinetic += 0.5*rho[]*dv ()*speed2;
  }
  return kinetic;
}

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
    EXTERNAL_FORCE_X,
    EXTERNAL_FORCE_Y,
    EXTERNAL_FORCE_Z,
    EXTERNAL_POWER,
    IBM_FORCE_X,
    IBM_FORCE_Y,
    IBM_FORCE_Z,
    IBM_POWER,
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
  // nforce is the current IBM reaction on the fibre.  gravity is the imposed
  // body-force density, and is the only force external to the coupled system.
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
    v[EXTERNAL_FORCE_X] += w*ibval (gravity.x);
    v[EXTERNAL_FORCE_Y] += w*ibval (gravity.y);
    v[EXTERNAL_FORCE_Z] += w*ibval (gravity.z);
    v[IBM_FORCE_X] += w*ibval (nforce.x);
    v[IBM_FORCE_Y] += w*ibval (nforce.y);
    v[IBM_FORCE_Z] += w*ibval (nforce.z);
    v[EXTERNAL_POWER] += w*(ibval (gravity.x)*ibval (nvel.x) +
                            ibval (gravity.y)*ibval (nvel.y) +
                            ibval (gravity.z)*ibval (nvel.z));
    v[IBM_POWER] += w*(ibval (nforce.x)*ibval (nvel.x) +
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
      fprintf (fp, "i,t,mesh_id,y_min,y_first_node,y_last_node,external_force_x,external_force_y,external_force_z,external_power,ibm_force_x,ibm_force_y,ibm_force_z,ibm_power\n");
    for (int m = 0; m < ibmm.nm; ++m) {
      double *v = values + STATS_STRIDE*m;
      double y_first = v[Y_FIRST_COUNT] > 0. ?
        v[Y_FIRST_SUM]/v[Y_FIRST_COUNT] : NAN;
      double y_last = v[Y_LAST_COUNT] > 0. ?
        v[Y_LAST_SUM]/v[Y_LAST_COUNT] : NAN;
      double y_min = minimum_y[m] < HUGE_VAL ? minimum_y[m] : NAN;
      fprintf (fp, "%d,%.17g,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
               iteration, time, m, y_min, y_first, y_last,
               v[EXTERNAL_FORCE_X], v[EXTERNAL_FORCE_Y], v[EXTERNAL_FORCE_Z],
               v[EXTERNAL_POWER], v[IBM_FORCE_X], v[IBM_FORCE_Y],
               v[IBM_FORCE_Z], v[IBM_POWER]);
    }
    fclose (fp);

    snprintf (filename, sizeof(filename), "%s/energy-diagnostics.csv", base_path);
    fp = fopen (filename, iteration == 0 ? "w" : "a");
    if (!fp) { perror (filename); exit (1); }
    if (iteration == 0)
      fprintf (fp, "i,t,dt,external_power,external_work,fluid_ibm_work,fibre_ibm_power,fibre_ibm_work,ibm_exchange_residual\n");
    bool valid_energy = isfinite (energy_step_dt) && energy_step_dt > 0. &&
      isfinite (external_power) && isfinite (fluid_ibm_work) &&
      isfinite (fibre_ibm_power);
    double external_work = valid_energy ? energy_step_dt*external_power : NAN;
    double fibre_ibm_work = valid_energy ? energy_step_dt*fibre_ibm_power : NAN;
    double residual = valid_energy ? fluid_ibm_work + fibre_ibm_work : NAN;
    fprintf (fp, "%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
             iteration, time, valid_energy ? energy_step_dt : NAN,
             valid_energy ? external_power : NAN, external_work,
             valid_energy ? fluid_ibm_work : NAN,
             valid_energy ? fibre_ibm_power : NAN, fibre_ibm_work, residual);
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

/** Snapshot fluid energy after alpha viscosity and immediately before MDF. */
event interface_force_richardson (i++, last) {
  fluid_kinetic_before_mdf = settling_fluid_kinetic_energy ();
}

/** Measure the MDF-only fluid energy change before beta viscosity runs.
 * nforce is the IBM reaction on the fibre produced by the inherited
 * interface_force_richardson event.  The equal-and-opposite Eulerian force is
 * already reflected in u, so the sum below is the interface exchange residual.
 */
event beta_viscous_term (i++, last) {
  double fluid_kinetic_after_mdf = settling_fluid_kinetic_energy ();
  fluid_ibm_work = fluid_kinetic_after_mdf - fluid_kinetic_before_mdf;
  double powers[2] = {0., 0.};
  foreach_ibnode_per_ibmesh () {
#if _MPI
    if (node->pid != pid ()) continue;
#endif
    double w = ibval (nweight);
    powers[0] += w*(ibval (nforce.x)*ibval (nvel.x) +
                     ibval (nforce.y)*ibval (nvel.y) +
                     ibval (nforce.z)*ibval (nvel.z));
    powers[1] += w*(ibval (gravity.x)*ibval (nvel.x) +
                     ibval (gravity.y)*ibval (nvel.y) +
                     ibval (gravity.z)*ibval (nvel.z));
  }
#if _MPI
  MPI_Allreduce (MPI_IN_PLACE, powers, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
  fibre_ibm_power = powers[0];
  external_power = powers[1];
  energy_step_dt = dt;
}

/** Report only after the full solver step; a standalone i-event declared by
 * the template would otherwise be scheduled immediately after advance_interface.
 */
event end_timestep (i++, last) {
  if (i % stats_interval == 0 && t <= generic_case.time_end)
    settling_fibre_stats (i, t);
}
