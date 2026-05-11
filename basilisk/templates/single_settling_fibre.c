#include "library/ibm/IBMeshManager.h"
#include "library/ibm/navier-stokes/centered-mdf.h"
#include "library/elff/elff.h"
#include "library/ibm/IBOutput.h"
#include "library/io/output-vtk.h"
#include "library/io/params/params-cli.h"

/* Default simulations parameters */

double dt_fluid = 0.0005;
double L_fluid = 16.;
double U0 = 0.;

int maxlevel = 10;
int minlevel = 4;
int ibmlevel = 11;

/*
 * Marchetti-style single-fibre control parameters.
 *
 * The length scale is the total fibre length L = b_length. We use
 * rho_f = 1 and g = 1 as the dimensional reference values for this
 * Basilisk setup, then derive the beam and fluid parameters that produce
 * the requested dimensionless groups.
 *
 *   marchetti_aspect_ratio = L/d = ell/a
 *   marchetti_Be           = Delta rho A_f g L^3 / EI
 *   marchetti_density_ratio = Delta rho / rho_f
 *   marchetti_Ga           = sqrt((Delta rho/rho_f) g L^3) / nu
 */
double marchetti_aspect_ratio = 30.;
double marchetti_Be = 200.;
double marchetti_density_ratio = 0.1;
double marchetti_Ga = 40.;

double b_length = 1.;
int b_nodes = 65;
double b_r = 1e3;
double b_theta = 0.01;

/* Diagnostic experiment controls */
double experiment_mass_scale = 1.;
double experiment_gravity_scale = 1.;
double experiment_gravity_ramp_time = 0.0;
double experiment_ib_force_relaxation = 0.4;
int experiment_ib_richardson_iters = 3;
double experiment_t_end = 500.0;
int experiment_stats_interval = 1;
double experiment_output_interval = 0.05;

char *base_path = "single_settling_fibre_output";

/* Derived parameters */

#define marchetti_rho_f (1.)
#define marchetti_g (1.)

#define b_diameter (b_length / marchetti_aspect_ratio)
#define b_radius (0.5 * b_diameter)
#define b_area (pi * b_radius * b_radius)
#define b_second_moment (pi * b_radius * b_radius * b_radius * b_radius / 4.)

#define b_rho_f (marchetti_rho_f)
#define b_delta_rho (marchetti_density_ratio * b_rho_f)
#define b_rho_s (b_rho_f + b_delta_rho)
#define b_mu (experiment_mass_scale * b_rho_s * b_area)
#define b_submerged_weight_per_length (b_delta_rho * b_area * marchetti_g)
#define b_EI                                                                   \
  (b_submerged_weight_per_length * b_length * b_length * b_length /            \
   marchetti_Be)
#define b_young_modulus (b_EI / b_second_moment)
#define b_gravity (b_submerged_weight_per_length)
#define fluid_nu                                                               \
  (sqrt(marchetti_density_ratio * marchetti_g * b_length * b_length *          \
        b_length) /                                                            \
   marchetti_Ga)
#define fluid_dynamic_viscosity (b_rho_f * fluid_nu)
#define marchetti_Re                                                           \
  (sqrt(marchetti_density_ratio * marchetti_g * b_length * b_length *          \
        b_length) *                                                            \
   b_diameter / fluid_nu)

/* Additional fields */

face vector muv[];

static void write_marchetti_parameters(FILE *fp) {
  fprintf(fp, "[marchetti]\n");
  fprintf(fp, "aspect_ratio_L_over_d = %.17g\n", marchetti_aspect_ratio);
  fprintf(fp, "Be = %.17g\n", marchetti_Be);
  fprintf(fp, "density_ratio_Delta_rho_over_rho_f = %.17g\n",
          marchetti_density_ratio);
  fprintf(fp, "Ga = %.17g\n", marchetti_Ga);

  fprintf(fp, "\n[derived]\n");
  fprintf(fp, "L = %.17g\n", b_length);
  fprintf(fp, "d = %.17g\n", b_diameter);
  fprintf(fp, "a = %.17g\n", b_radius);
  fprintf(fp, "A_f = %.17g\n", b_area);
  fprintf(fp, "I = %.17g\n", b_second_moment);
  fprintf(fp, "rho_f = %.17g\n", b_rho_f);
  fprintf(fp, "rho_s = %.17g\n", b_rho_s);
  fprintf(fp, "Delta_rho = %.17g\n", b_delta_rho);
  fprintf(fp, "mu_s = %.17g\n", b_mu);
  fprintf(fp, "W = %.17g\n", b_submerged_weight_per_length);
  fprintf(fp, "EI = %.17g\n", b_EI);
  fprintf(fp, "E = %.17g\n", b_young_modulus);
  fprintf(fp, "nu = %.17g\n", fluid_nu);
  fprintf(fp, "mu_f = %.17g\n", fluid_dynamic_viscosity);
  fprintf(fp, "Re_d = %.17g\n", marchetti_Re);
  fprintf(fp, "b_submerged_weight_per_length = %.17g\n",
          b_submerged_weight_per_length);

  fprintf(fp, "\n[experiment]\n");
  fprintf(fp, "mass_scale = %.17g\n", experiment_mass_scale);
  fprintf(fp, "gravity_scale = %.17g\n", experiment_gravity_scale);
  fprintf(fp, "gravity_ramp_time = %.17g\n", experiment_gravity_ramp_time);
  fprintf(fp, "ib_force_relaxation = %.17g\n", experiment_ib_force_relaxation);
  fprintf(fp, "ib_richardson_iters = %d\n", experiment_ib_richardson_iters);
  fprintf(fp, "t_end = %.17g\n", experiment_t_end);
  fprintf(fp, "stats_interval = %d\n", experiment_stats_interval);
  fprintf(fp, "output_interval = %.17g\n", experiment_output_interval);
}

/* Quiescent tank boundary conditions */

u.n[left] = dirichlet(0.);
u.t[left] = dirichlet(0.);
p[left] = neumann(0.);
pf[left] = neumann(0.);

u.n[right] = dirichlet(0.);
u.t[right] = dirichlet(0.);
p[right] = neumann(0.);
pf[right] = neumann(0.);

u.n[top] = dirichlet(0.);
u.t[top] = dirichlet(0.);
p[top] = neumann(0.);
pf[top] = neumann(0.);

u.n[bottom] = dirichlet(0.);
u.t[bottom] = dirichlet(0.);
p[bottom] = neumann(0.);
pf[bottom] = neumann(0.);

int main(int argc, char **argv) {
  /* Here we register runtime options for the simulation */
  input_file_register_option("basilisk.fluid", L_fluid, PARAM_VALUE_DOUBLE);
  input_file_register_option("basilisk.fluid", dt_fluid, PARAM_VALUE_DOUBLE);
  input_file_register_option("basilisk.fluid", minlevel, PARAM_VALUE_INT);
  input_file_register_option("basilisk.fluid", maxlevel, PARAM_VALUE_INT);
  input_file_register_option("basilisk.fluid", ibmlevel, PARAM_VALUE_INT);
  input_file_register_option("basilisk.output", base_path, PARAM_VALUE_STRING);
  input_file_register_option("marchetti", marchetti_aspect_ratio,
                             PARAM_VALUE_DOUBLE);
  input_file_register_option("marchetti", marchetti_Be, PARAM_VALUE_DOUBLE);
  input_file_register_option("marchetti", marchetti_density_ratio,
                             PARAM_VALUE_DOUBLE);
  input_file_register_option("marchetti", marchetti_Ga, PARAM_VALUE_DOUBLE);
  input_file_register_option("beam", b_length, PARAM_VALUE_DOUBLE);
  input_file_register_option("beam", b_nodes, PARAM_VALUE_INT);
  input_file_register_option("beam", b_r, PARAM_VALUE_DOUBLE);
  input_file_register_option("beam", b_theta, PARAM_VALUE_DOUBLE);
  input_file_register_option("experiment", experiment_mass_scale,
                             PARAM_VALUE_DOUBLE);
  input_file_register_option("experiment", experiment_gravity_scale,
                             PARAM_VALUE_DOUBLE);
  input_file_register_option("experiment", experiment_gravity_ramp_time,
                             PARAM_VALUE_DOUBLE);
  input_file_register_option("experiment", experiment_ib_force_relaxation,
                             PARAM_VALUE_DOUBLE);
  input_file_register_option("experiment", experiment_ib_richardson_iters,
                             PARAM_VALUE_INT);
  input_file_register_option("experiment", experiment_t_end,
                             PARAM_VALUE_DOUBLE);
  input_file_register_option("experiment", experiment_stats_interval,
                             PARAM_VALUE_INT);
  input_file_register_option("experiment", experiment_output_interval,
                             PARAM_VALUE_DOUBLE);

  /* Here we parse the options given through the command-line */
  int input_file_parse_result = input_file_parse_cli(argc, argv);
  if (input_file_parse_result != 0)
    return input_file_parse_result > 0 ? 0 : input_file_parse_result;
  else {
    input_file_print_options();
    input_file_apply_options();
  }

  /* Setting relevant parameters for basilisk */
  L0 = L_fluid;
  origin(-L_fluid / 2., -L_fluid / 2., -L_fluid / 2.);
#if TREE
  N = 1 << minlevel;
#else
  N = 1 << ibmlevel;
#endif
  DT = dt_fluid;
  mu = muv;
  ib_force_relaxation = experiment_ib_force_relaxation;
  ib_richardson_iters = experiment_ib_richardson_iters;
  display_control(marchetti_Re, 10, 1000);

  if (pid() == 0) {
    write_marchetti_parameters(stderr);

    create_path(base_path);
    char fname[4096];
    snprintf(fname, sizeof(fname), "%s/marchetti_parameters.txt", base_path);
    FILE *fp = fopen(fname, "w");
    if (fp) {
      write_marchetti_parameters(fp);
      fclose(fp);
    } else {
      fprintf(stderr, "warning: failed to open %s for writing\n", fname);
    }
  }

  run();
}

/* Reynolds control */
event properties(i++) {
  foreach_face() muv.x[] = fm.x[] * fluid_dynamic_viscosity;
}

/* Beam setup */
event init(i = 0) {
  int m_id = ibmeshmanager_add_mesh();
  // ib_euler_beam_bcs_t b_bcs = elff_euler_beam_bcs_theta_pin((coord){0});
  ib_euler_beam_bcs_t b_bcs = {
      .end = {IB_EULER_BEAM_BC_LEFT, IB_EULER_BEAM_BC_RIGHT},
      .type = {IB_EULER_BEAM_BC_FREE, IB_EULER_BEAM_BC_FREE},
      .vals = {0},
  };

  IBMeshModel beam_model = elff_euler_beam_addm_new_theta(
      b_length, b_EI, b_mu, b_nodes, b_r, b_theta, b_bcs);
  // IBMeshModel beam_model = elff_euler_beam_addm_new_theta(b_length,
  // 1e-6, 1.1, b_nodes, b_r, b_theta, b_bcs);

  ibmeshmanager_set_model(m_id, beam_model);

  foreach_ibnode_per_ibmesh() {
    mesh->depth = ibmlevel;
    node->depth = ibmlevel;
    ibval(gravity.x) = 0.;
    ibval(gravity.y) = 0.;
  }

  if (!restore_handler(base_path)) {
    foreach ()
      foreach_dimension() u.x[] = 0.;
  } else {
  }
}

event body_force(i++) {
  double ramp = 1.;
  if (experiment_gravity_ramp_time > 0.)
    ramp = min(t / experiment_gravity_ramp_time, 1.);

  foreach_ibnode_per_ibmesh() {
    ibval(gravity.x) = 0.;
    ibval(gravity.y) =
        -experiment_gravity_scale * ramp * b_submerged_weight_per_length;
  }
}

event logfile(i++) {
  if (pid() == 0)
    fprintf(stderr, "[info] %d %g\n", i, t);
}

event statsfile(i += experiment_stats_interval; t <= experiment_t_end) {
  double x_first = 0.;
  double y_first = 0.;
  double x_last = 0.;
  double y_last = 0.;
  double x_com = 0.;
  double y_com = 0.;
  double w_sum = 0.;
  int node_count = 0;

  foreach_ibnode() {
    const double w = ibval(nweight);
    if (node_count == 0) {
      x_first = ibval(npos.x);
      y_first = ibval(npos.y);
    }

    x_last = ibval(npos.x);
    y_last = ibval(npos.y);
    x_com += w * ibval(npos.x);
    y_com += w * ibval(npos.y);
    w_sum += w;
    node_count++;
  }

  if (w_sum > 0.) {
    x_com /= w_sum;
    y_com /= w_sum;
  }

  const double end_to_end = sqrt(sq(x_last - x_first) + sq(y_last - y_first));

  if (pid() == 0) {
    FILE *fp = NULL;
    if (!base_path) {
      fprintf(stderr, "warning: base_path is NULL; skipping tip output\n");
      return 0;
    }
    create_path(base_path);
    char fname[4096];
    snprintf(fname, sizeof(fname), "%s/tip.txt", base_path);
    fp = fopen(fname, i == 0 ? "w" : "a");
    if (!fp) {
      fprintf(stderr, "warning: failed to open %s for append\n", fname);
      return 0;
    }
    if (i == 0)
      fprintf(fp,
              "# i t x_first y_first x_last y_last x_com y_com end_to_end\n");
    fprintf(fp, "%d %g %g %g %g %g %g %g %g\n", i, t, x_first, y_first, x_last,
            y_last, x_com, y_com, end_to_end);
    fclose(fp);
  }
}

event output(i += 500; t <= experiment_t_end) {
  scalar omega[];
  vorticity(u, omega);
#if TREE
  output_hdf_htg({omega, p}, {u, ibmf}, base_path);
#else
  output_hdf_imagedata({omega, p}, {u, ibmf}, base_path);
#endif
  output_hdf_pd(NULL, (IBvector[]){eulvel, nforce, nvel, {{-1}}}, base_path);
}

#if TREE
event adapt(i++) {
  adapt_wavelet_ibm({u}, (double[]){3e-3, 3e-3, 3e-3}, maxlevel, minlevel);
}
#endif

event checkpoint_event(i++, last) {
  return checkpoint_handler(t, i, base_path);
}
