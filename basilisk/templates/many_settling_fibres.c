#include <math.h>

#ifndef IBMESHMANAGER_POOL_SIZE_BYTES
#define IBMESHMANAGER_POOL_SIZE_BYTES (1 << 24)
#endif

#include "library/ibm/IBMeshManager.h"
#include "library/ibm/navier-stokes/centered-mdf.h"
#include "library/elff/elff.h"
#include "library/ibm/IBOutput.h"
#include "library/io/output-vtk.h"
#include "library/io/params/params-cli.h"

#include "lambda2.h"

#if dimension != 3
#error "many_settling_fibres_huang is configured for grid/octree.h"
#endif

/* Default simulation parameters */

double dt_fluid = 0.001;
double L_fluid = 32.;
double U0 = 0.;

int maxlevel = 12;
int minlevel = 4;
int ibmlevel = 12;

/*
 * Banaei et al. (2020) nondimensional groups:
 *   rp    = L/d
 *   gamma = EI/(rho_l r g L^3), with rho_l = Delta rho A_f
 *   r     = Delta rho/rho_0
 *   Ga    = sqrt(r g L^3)/nu
 */

double banaei_rp = 30.;
double banaei_gamma = 0.01;
double banaei_r = 0.1;
double banaei_Ga = 40.;

double b_length = 1.;
int b_nodes = 67;
double b_theta = 0.00;
int b_pid = 0;

/* Fibre cloud controls */
int fibre_count = 1000;
int fibre_seed = 20260727;
double fibre_cloud_width = 0.;
double fibre_cloud_height = 0.;
double fibre_center_y = 0.;
double fibre_margin = 0.;
double fibre_min_center_spacing_ratio = 1.25;
int fibre_random_orientation = 1;
int fibre_random_azimuth = 1;

/* Nonphysical experiment controls */
double experiment_ib_force_relaxation = 0.8;
int experiment_ib_richardson_iters = 2;
int experiment_huang_implicit_bending = 1;
double experiment_adapt_velocity_tolerance = 1e-3;

int experiment_stats_interval = 1;
double experiment_output_interval = 1.0;
double experiment_output_iters = 1;

double experiment_t_end = 500.0;

char *base_path = "many_settling_fibres_huang_output";

/* Derived parameters */
#define b_rho_0 (1.)
#define b_g (1.)

#define b_diameter (b_length / banaei_rp)
#define b_area (pi * b_diameter * b_diameter / 4.)
#define b_linear_density_difference (banaei_r * b_rho_0 * b_area)
#define b_mu ((1. + banaei_r) * b_rho_0 * b_area)
#define b_submerged_weight_per_length (b_linear_density_difference * b_g)
#define b_EI                                                                 \
  (banaei_gamma * banaei_r * b_submerged_weight_per_length * b_length *       \
   b_length * b_length)
#define b_gravity -(b_submerged_weight_per_length)
#define b_ds (b_length / (b_nodes - 1))

#define fluid_nu                                                              \
  (sqrt(banaei_r * b_g * b_length * b_length * b_length) / banaei_Ga)
#define fluid_dynamic_viscosity (b_rho_0 * fluid_nu)
#define fluid_velocity_scale (sqrt(banaei_r * b_g * b_length))
#define h_ibm (L_fluid / (1 << ibmlevel))

/* Additional fields */

face vector muv[];

static unsigned rng_state = 1u;

static double uniform_random(void) {
  rng_state = 1664525u * rng_state + 1013904223u;
  return (double)(rng_state >> 8) / 16777216.;
}

static double uniform_range(double lo, double hi) {
  return lo + (hi - lo) * uniform_random();
}

static double clamp_double(double value, double lo, double hi) {
  if (hi < lo)
    return lo;
  return min(max(value, lo), hi);
}

static double effective_fibre_margin(void) {
  const double requested_margin =
    fibre_margin > 0. ? fibre_margin : 0.5 * b_length + 2. * h_ibm;
  return min(requested_margin, 0.45 * L_fluid);
}

static coord random_fibre_direction(void) {
  if (fibre_random_orientation) {
    const double azimuth = uniform_range(0., 2. * pi);
    const double direction_z = uniform_range(-1., 1.);
    const double direction_xy = sqrt(max(0., 1. - sq(direction_z)));
    coord direction = {
      direction_xy * cos(azimuth),
      direction_xy * sin(azimuth),
      direction_z
    };
    return direction;
  }

  const double azimuth =
    fibre_random_azimuth ? uniform_range(0., 2. * pi) : 0.;
  const double horizontal = cos(b_theta);

  coord direction = {
    horizontal * cos(azimuth),
    sin(b_theta),
    horizontal * sin(azimuth)
  };

  const double norm =
    sqrt(sq(direction.x) + sq(direction.y) + sq(direction.z));
  foreach_dimension()
    direction.x /= norm;

  return direction;
}

static coord random_fibre_center(double margin) {
  const double half_allowed = max(0., 0.5 * L_fluid - margin);
  const double half_width =
    fibre_cloud_width > 0. ? min(0.5 * fibre_cloud_width, half_allowed)
                           : half_allowed;

  const double y_center =
    clamp_double(fibre_center_y, -half_allowed, half_allowed);
  const double half_height_allowed =
    min(y_center + half_allowed, half_allowed - y_center);
  const double half_height =
    fibre_cloud_height > 0.
      ? min(0.5 * fibre_cloud_height, max(0., half_height_allowed))
      : half_allowed;
  const double y_position =
    fibre_cloud_height > 0. ? uniform_range(y_center - half_height,
                                            y_center + half_height)
                            : uniform_range(-half_allowed, half_allowed);

  coord center = {
    uniform_range(-half_width, half_width),
    y_position,
    uniform_range(-half_width, half_width)
  };
  return center;
}

static int center_is_separated(coord center,
                               coord *centers,
                               int placed,
                               double min_center_spacing) {
  if (min_center_spacing <= 0.)
    return 1;

  for (int n = 0; n < placed; ++n) {
    const double distance2 =
      sq(center.x - centers[n].x) + sq(center.y - centers[n].y) +
      sq(center.z - centers[n].z);
    if (distance2 < sq(min_center_spacing))
      return 0;
  }

  return 1;
}

static void add_settling_fibres(void) {
  ib_euler_beam_bcs_t b_bcs;
  const double margin = effective_fibre_margin();
  const double min_center_spacing =
    fibre_min_center_spacing_ratio * b_length;
  coord centers[fibre_count];

  b_bcs.end[0] = IB_EULER_BEAM_BC_LEFT;
  b_bcs.end[1] = IB_EULER_BEAM_BC_RIGHT;
  b_bcs.type[0] = IB_EULER_BEAM_BC_FREE;
  b_bcs.type[1] = IB_EULER_BEAM_BC_FREE;

  rng_state = (unsigned)(fibre_seed ? fibre_seed : 1);

  for (int n = 0; n < fibre_count; ++n) {
    coord center = {0.};
    coord direction = {1., 0., 0.};
    int accepted = 0;

    for (int attempt = 0; attempt < 10000; ++attempt) {
      center = random_fibre_center(margin);
      direction = random_fibre_direction();
      if (center_is_separated(center, centers, n, min_center_spacing)) {
        accepted = 1;
        break;
      }
    }

    if (!accepted && pid() == 0)
      fprintf(stderr,
              "warning: accepted fibre %d without satisfying centre spacing\n",
              n);

    centers[n] = center;

    coord start = {
      center.x - 0.5 * b_length * direction.x,
      center.y - 0.5 * b_length * direction.y,
      center.z - 0.5 * b_length * direction.z
    };

    int mesh_id = ibmeshmanager_add_mesh();
    IBMeshModel beam_model =
      elff_euler_beam_huang_new_direction(
        b_length, b_EI, b_mu, b_nodes, direction, b_bcs, start, b_pid);
    ib_euler_beam_huang_set_implicit_bending(
      (ib_euler_beam_huang_t)beam_model.ctx,
      experiment_huang_implicit_bending);

    vertex_t initial_velocity[b_nodes];
    for (int ni = 0; ni < b_nodes; ++ni)
      initial_velocity[ni] = (vertex_t){0., 0., 0.};
    ib_euler_beam_huang_set_initial_velocity(
      (ib_euler_beam_huang_t)beam_model.ctx,
      initial_velocity, b_nodes);

    ibmeshmanager_set_model(mesh_id, beam_model);
  }

  foreach_ibnode_per_ibmesh() {
    mesh->depth = ibmlevel;
    node->depth = ibmlevel;
    foreach_dimension()
      ibval(gravity.x) = 0.;
    ibval(gravity.y) = b_gravity;
  }
}

static void write_settling_stats(int iter,
                                 double time,
                                 double eulerian_ibm_force_y) {
  const int meshes = ibmm.nm;
  if (meshes <= 0)
    return;

  double measure[meshes];
  double sum_x[meshes];
  double sum_y[meshes];
  double sum_z[meshes];
  double min_y[meshes];
  double max_y[meshes];
  double sum_velocity_y[meshes];
  double nodal_force_y = 0.;
  double measure_sum = 0.;

  for (int m = 0; m < meshes; ++m) {
    measure[m] = 0.;
    sum_x[m] = 0.;
    sum_y[m] = 0.;
    sum_z[m] = 0.;
    min_y[m] = Y0 + L0;
    max_y[m] = Y0;
    sum_velocity_y[m] = 0.;
  }

  foreach_ibnode_per_ibmesh() {
    const int m = (int)mesh_id;
    const double w = ibval(nweight);

    measure[m] += w;
    sum_x[m] += ibval(npos.x) * w;
    sum_y[m] += ibval(npos.y) * w;
    sum_z[m] += ibval(npos.z) * w;
    sum_velocity_y[m] += ibval(nvel.y) * w;
    min_y[m] = min(min_y[m], ibval(npos.y));
    max_y[m] = max(max_y[m], ibval(npos.y));
    nodal_force_y += ibval(nforce.y) * w;
    measure_sum += w;
  }

  double mean_x = 0.;
  double mean_y = 0.;
  double mean_z = 0.;
  double mean_velocity_y = 0.;
  double global_min_y = Y0 + L0;
  double global_max_y = Y0;

  for (int m = 0; m < meshes; ++m) {
    mean_x += sum_x[m];
    mean_y += sum_y[m];
    mean_z += sum_z[m];
    mean_velocity_y += sum_velocity_y[m];
    global_min_y = min(global_min_y, min_y[m]);
    global_max_y = max(global_max_y, max_y[m]);
  }

  if (measure_sum > 0.) {
    mean_x /= measure_sum;
    mean_y /= measure_sum;
    mean_z /= measure_sum;
    mean_velocity_y /= measure_sum;
  }

  const double hydro_force_y = -eulerian_ibm_force_y;
  const double gravity_force_y = b_gravity * measure_sum;
  const double force_balance_y = hydro_force_y + gravity_force_y;

  if (pid() != 0)
    return;

  fprintf(stderr,
          "diag i=%d t=%g fibres=%d nodes=%zu mean_y=%g min_y=%g "
          "mean_vy=%g hydro_y=%g gravity_y=%g balance_y=%g\n",
          iter, time, meshes, ibmm.pool.active.size, mean_y, global_min_y,
          mean_velocity_y, hydro_force_y, gravity_force_y, force_balance_y);

  if (!base_path) {
    fprintf(stderr, "warning: base_path is NULL; skipping settling stats\n");
    return;
  }

  if (create_path(base_path) != 0) {
    fprintf(stderr,
            "warning: failed to create output path %s; skipping stats\n",
            base_path);
    return;
  }

  char fname[4096];
  snprintf(fname, sizeof(fname), "%s/settling-stats.csv", base_path);
  FILE *fp = fopen(fname, iter == 0 ? "w" : "a");
  if (!fp) {
    fprintf(stderr, "warning: failed to open %s for append\n", fname);
    return;
  }

  if (iter == 0)
    fprintf(fp,
            "i,t,banaei_Ga,banaei_gamma,banaei_r,banaei_rp,B,"
            "fibres,nodes_per_fibre,total_nodes,mean_x,mean_y,mean_z,"
            "min_y,max_y,mean_velocity_y,measure_sum,hydro_force_y,"
            "eulerian_ibm_force_y,nodal_force_y,gravity_force_y,"
            "force_balance_y\n");

  fprintf(fp,
          "%d,%g,%g,%g,%g,%g,%g,%d,%d,%zu,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g\n",
          iter, time, banaei_Ga, banaei_gamma, banaei_r, banaei_rp,
          1. / (banaei_gamma * banaei_r), meshes, b_nodes,
          ibmm.pool.active.size, mean_x, mean_y, mean_z, global_min_y,
          global_max_y, mean_velocity_y, measure_sum, hydro_force_y,
          eulerian_ibm_force_y, nodal_force_y, gravity_force_y,
          force_balance_y);
  fclose(fp);

  snprintf(fname, sizeof(fname), "%s/fibre-centres.csv", base_path);
  fp = fopen(fname, iter == 0 ? "w" : "a");
  if (!fp) {
    fprintf(stderr, "warning: failed to open %s for append\n", fname);
    return;
  }

  if (iter == 0)
    fprintf(fp,
            "i,t,fibre_id,center_x,center_y,center_z,min_y,max_y,"
            "mean_velocity_y,measure\n");

  for (int m = 0; m < meshes; ++m) {
    const double inv_measure = measure[m] > 0. ? 1. / measure[m] : 0.;
    fprintf(fp,
            "%d,%g,%d,%g,%g,%g,%g,%g,%g,%g\n",
            iter, time, m, sum_x[m] * inv_measure, sum_y[m] * inv_measure,
            sum_z[m] * inv_measure, min_y[m], max_y[m],
            sum_velocity_y[m] * inv_measure, measure[m]);
  }
  fclose(fp);
}

/* Triply-periodic box; periodic faces are set after the domain origin. */

int main(int argc, char **argv) {
  input_file_register_option("basilisk.fluid", L_fluid, PARAM_VALUE_DOUBLE);
  input_file_register_option("basilisk.fluid", dt_fluid, PARAM_VALUE_DOUBLE);
  input_file_register_option("basilisk.fluid", minlevel, PARAM_VALUE_INT);
  input_file_register_option("basilisk.fluid", maxlevel, PARAM_VALUE_INT);
  input_file_register_option("basilisk.fluid", ibmlevel, PARAM_VALUE_INT);
  input_file_register_option("basilisk.output", base_path, PARAM_VALUE_STRING);
  input_file_register_option_named("banaei", "rp", banaei_rp,
                                   PARAM_VALUE_DOUBLE);
  input_file_register_option_named("banaei", "gamma", banaei_gamma,
                                   PARAM_VALUE_DOUBLE);
  input_file_register_option_named("banaei", "r", banaei_r,
                                   PARAM_VALUE_DOUBLE);
  input_file_register_option_named("banaei", "Ga", banaei_Ga,
                                   PARAM_VALUE_DOUBLE);
  input_file_register_option("beam", b_length, PARAM_VALUE_DOUBLE);
  input_file_register_option("beam", b_nodes, PARAM_VALUE_INT);
  input_file_register_option("beam", b_theta, PARAM_VALUE_DOUBLE);
  input_file_register_option("fibre", fibre_count, PARAM_VALUE_INT);
  input_file_register_option("fibre", fibre_seed, PARAM_VALUE_INT);
  input_file_register_option("fibre", fibre_cloud_width, PARAM_VALUE_DOUBLE);
  input_file_register_option("fibre", fibre_cloud_height, PARAM_VALUE_DOUBLE);
  input_file_register_option("fibre", fibre_center_y, PARAM_VALUE_DOUBLE);
  input_file_register_option("fibre", fibre_margin, PARAM_VALUE_DOUBLE);
  input_file_register_option("fibre", fibre_min_center_spacing_ratio,
                             PARAM_VALUE_DOUBLE);
  input_file_register_option("fibre", fibre_random_orientation,
                             PARAM_VALUE_INT);
  input_file_register_option("fibre", fibre_random_azimuth, PARAM_VALUE_INT);
  input_file_register_option("experiment", experiment_ib_force_relaxation,
                             PARAM_VALUE_DOUBLE);
  input_file_register_option("experiment", experiment_ib_richardson_iters,
                             PARAM_VALUE_INT);
  input_file_register_option("experiment", experiment_huang_implicit_bending,
                             PARAM_VALUE_INT);
  input_file_register_option("experiment",
                             experiment_adapt_velocity_tolerance,
                             PARAM_VALUE_DOUBLE);
  input_file_register_option("experiment", experiment_t_end,
                             PARAM_VALUE_DOUBLE);
  input_file_register_option("experiment", experiment_stats_interval,
                             PARAM_VALUE_INT);
  input_file_register_option("experiment", experiment_output_interval,
                             PARAM_VALUE_DOUBLE);

  int input_file_parse_result = input_file_parse_cli(argc, argv);
  if (input_file_parse_result != 0)
    return input_file_parse_result > 0 ? 0 : input_file_parse_result;

  input_file_print_options();
  input_file_apply_options();

  if (L_fluid <= 0. || dt_fluid <= 0. || experiment_t_end <= 0. ||
      experiment_output_interval <= 0.) {
    fprintf(stderr, "fluid scales and time controls must be positive\n");
    return 2;
  }
  if (banaei_rp <= 0. || banaei_gamma <= 0. || banaei_r <= 0. ||
      banaei_Ga <= 0. || b_length <= 0. || b_nodes <= 5) {
    fprintf(stderr, "beam geometry and Banaei groups must be positive\n");
    return 2;
  }
  if (fibre_count <= 0 || fibre_cloud_width < 0. ||
      fibre_cloud_height < 0. || fibre_min_center_spacing_ratio < 0.) {
    fprintf(stderr, "fibre count and cloud controls are invalid\n");
    return 2;
  }
  if (minlevel < 0 || maxlevel < minlevel || ibmlevel < minlevel ||
      maxlevel < ibmlevel) {
    fprintf(stderr, "levels must satisfy 0 <= minlevel <= ibmlevel <= maxlevel\n");
    return 2;
  }
  if (experiment_stats_interval <= 0 ||
      experiment_adapt_velocity_tolerance <= 0.) {
    fprintf(stderr, "experiment intervals and tolerances must be positive\n");
    return 2;
  }

  L0 = L_fluid;
  origin(-L_fluid / 2., -L_fluid / 2., -L_fluid / 2.);
  foreach_dimension()
    periodic(right);
#if TREE
  N = 1 << minlevel;
#else
  N = 1 << ibmlevel;
#endif
  DT = dt_fluid;
  mu = muv;
  ib_force_relaxation = experiment_ib_force_relaxation;
  ib_richardson_iters = experiment_ib_richardson_iters;
  display_control(banaei_Ga, 1, 1000);

  if (pid() == 0)
    fprintf(stderr,
            "many-fibre settling: L0/L=%g Ga=%g gamma=%g r=%g rp=%g "
            "B=%g fibres=%d nodes=%d ds/dx=%g cloud_width/L=%g "
            "cloud_height/L=%g margin/L=%g min_center_spacing/L=%g "
            "random_orientation=%d nu=%g dt=%g t_end=%g\n",
            L_fluid / b_length, banaei_Ga, banaei_gamma, banaei_r,
            banaei_rp, 1. / (banaei_gamma * banaei_r), fibre_count,
            b_nodes, b_ds / h_ibm, fibre_cloud_width / b_length,
            fibre_cloud_height / b_length, effective_fibre_margin() / b_length,
            fibre_min_center_spacing_ratio, fibre_random_orientation, fluid_nu, dt_fluid,
            experiment_t_end);

  run();
}

event logfile(i += experiment_stats_interval; t <= experiment_t_end) {
  double eulerian_ibm_force_y = 0.;

  foreach(reduction(+:eulerian_ibm_force_y))
    eulerian_ibm_force_y += ibmf.y[] * dv();

  write_settling_stats(i, t, eulerian_ibm_force_y);
}

event properties(i++) {
  foreach_face()
    muv.x[] = fm.x[] * fluid_dynamic_viscosity;
}

event init(i = 0) {
  add_settling_fibres();

  if (!restore_handler(base_path)) {
    foreach()
      foreach_dimension()
        u.x[] = 0.;
#if TREE
    adapt_wavelet_ibm(NULL, NULL, 0, 1, all, true);
#endif
  }
}

event output(i += experiment_output_iters; t <= experiment_t_end) {
  scalar l2[], omega_z[];
  lambda2(u, l2);
  vorticity(u, omega_z);

#if TREE
  output_hdf_htg({l2, omega_z, p}, {u, ibmf}, base_path);
#else
  output_hdf_imagedata({l2, omega_z, p}, {u, ibmf}, base_path);
#endif
  output_hdf_pd(NULL, (IBvector[]){eulvel, nforce, nvel, {{-1}}}, base_path);
}

#if TREE
event adapt(i++) {
  const double val =
    experiment_adapt_velocity_tolerance * fluid_velocity_scale;
  adapt_wavelet_ibm({u}, (double[]){val, val, val}, maxlevel, minlevel);
}
#endif

event checkpoint_event(i++, last) {
  return checkpoint_handler(t, i, base_path);
}

event stop(t = experiment_t_end) {
  if (pid() == 0)
    fprintf(stderr, "finished many-fibre settling at t=%g i=%d\n", t, i);
  return 1;
}
