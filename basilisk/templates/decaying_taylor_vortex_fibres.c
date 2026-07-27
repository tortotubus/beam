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
#error "decaying_taylor_vortex_fibres is configured for grid/octree.h"
#endif

/*
 * Fibre nondimensionalization:
 *   L_c   = L_f              fibre length
 *   U_c   = U0               initial velocity amplitude
 *   T_c   = L_f/U0           fibre advection time
 *   rho_c = rho_f            fluid density
 *
 * For a circular fibre with material density ratio R = rho_s/rho_f and
 * aspect ratio AR = L_f/D,
 *
 *   q   = rho_b/(rho_f L_f^2) = R*pi/(4*AR^2),
 *   K_B = EI/(rho_b U_c^2 L_f^2).
 *
 * Suspension concentration is prescribed by
 *
 *   c = n*L_f^3 = N*L_f^3/V,
 *
 * where n is number density and V is the periodic box volume. Neglecting
 * fibre end caps, the corresponding solid volume fraction is
 * phi = c*pi/(4*AR^2).
 *
 * For fixed N and c, the cubic box side is
 *
 *   L_box = L_f*(N/c)^(1/3).
 *
 * The Taylor-Green length and time are L_TG = 1/k = L_box/(2*pi) and
 * T_TG = L_TG/U0. Reynolds, dt_fluid, t_end, and dt_output are specified
 * on these Taylor-Green scales and converted to simulation units below.
 */
char *base_path = "decaying_taylor_vortex_fibres_output";

double U0 = 1.;
double Reynolds = 200.;
double dt_fluid = 1e-2;
double t_end = 10.;
double dt_output = 0.1;

int maxlevel = 10;
int minlevel = 6;
int ibmlevel = 10;

/* Fibre parameters */
int fibre_count = 100;
int fibre_seed = 20260726;
double fibre_length = 1.;
double fibre_aspect_ratio = 30.;
/* Gives L_box/L_f = 20*pi for the default count. */
double fibre_nL3 = 0.5; //1. / (8. * pi * pi * pi);
/* Matches ds/dx = (1/66)/(16/2048) in single_settling_fibre_huang.c. */
double fibre_marker_spacing_ratio = 128. / 66.;
double fibre_density_ratio = 1.2;
double fibre_dimensionless_bending_rigidity = 1e-3;
double fibre_margin = 0.;
int fibre_pid = 0;

/* Coupling and diagnostics */
double experiment_ib_force_relaxation = 1.;
int experiment_ib_richardson_iters = 2;
int experiment_huang_implicit_bending = 1;
double experiment_adapt_velocity_tolerance = 0.05;
int stats_interval = 10;

#define L_fluid (fibre_length * cbrt(fibre_count / fibre_nL3))
#define taylor_k (2. * pi / L_fluid)
#define taylor_length_scale (1. / taylor_k)
#define taylor_time_scale (taylor_length_scale / U0)
#define fluid_density (1.)
#define fluid_nu (U0 * taylor_length_scale / Reynolds)
#define fluid_dynamic_viscosity (fluid_density * fluid_nu)
#define simulation_dt (dt_fluid * taylor_time_scale)
#define simulation_end_time (t_end * taylor_time_scale)
#define simulation_output_interval (dt_output * taylor_time_scale)
#define h_ibm (L_fluid / (1 << ibmlevel))
#define fibre_diameter (fibre_length / fibre_aspect_ratio)
#define fibre_area (pi * sq(fibre_diameter) / 4.)
#define fibre_line_density_ratio \
  (fibre_density_ratio * fibre_area / sq(fibre_length))
#define fibre_linear_density \
  (fibre_line_density_ratio * fluid_density * sq(fibre_length))
#define fibre_bending_stiffness                                      \
  (fibre_dimensionless_bending_rigidity * fibre_linear_density *    \
   sq(U0 * fibre_length))
#define fibre_volume_fraction \
  (fibre_nL3 * pi / (4. * sq(fibre_aspect_ratio)))
#define fibre_nodes                                                      \
  ((int)floor(fibre_length / (fibre_marker_spacing_ratio * h_ibm) + 0.5) + 1)
#define fibre_ds (fibre_length / (fibre_nodes - 1))

face vector muv[];

static unsigned rng_state = 1u;

static double uniform_random(void) {
  rng_state = 1664525u * rng_state + 1013904223u;
  return (double)(rng_state >> 8) / 16777216.;
}

static double uniform_range(double lo, double hi) {
  return lo + (hi - lo) * uniform_random();
}

static coord taylor_velocity(coord position) {
  const double kx = taylor_k * position.x;
  const double ky = taylor_k * position.y;
  const double kz = taylor_k * position.z;
  coord velocity = {
    -U0 * cos(kx) * sin(ky) * cos(kz),
     U0 * sin(kx) * cos(ky) * cos(kz),
     0.
  };
  return velocity;
}

static coord taylor_directional_velocity_gradient(coord position,
                                                   coord direction) {
  const double kx = taylor_k * position.x;
  const double ky = taylor_k * position.y;
  const double kz = taylor_k * position.z;
  const double scale = U0 * taylor_k;
  coord gradient = {
    scale * ( sin(kx) * sin(ky) * cos(kz) * direction.x
             -cos(kx) * cos(ky) * cos(kz) * direction.y
             +cos(kx) * sin(ky) * sin(kz) * direction.z),
    scale * ( cos(kx) * cos(ky) * cos(kz) * direction.x
             -sin(kx) * sin(ky) * cos(kz) * direction.y
             -sin(kx) * cos(ky) * sin(kz) * direction.z),
    0.
  };
  return gradient;
}

static void set_taylor_green_field(void) {
  foreach() {
    const double kx = taylor_k * x;
    const double ky = taylor_k * y;
    const double kz = taylor_k * z;
    u.x[] = -U0 * cos(kx) * sin(ky) * cos(kz);
    u.y[] =  U0 * sin(kx) * cos(ky) * cos(kz);
    u.z[] = 0.;
    p[] = -sq(U0) / 16. *
          (cos(2. * kx) + cos(2. * ky)) * (cos(2. * kz) + 2.);
  }
  boundary((scalar *){u, p});

  foreach()
    foreach_dimension()
      g.x[] = -(p[1] - p[-1]) / (2. * Delta);
  boundary((scalar *){g});
}

static void add_random_fibres(void) {
  double requested_margin;
  double margin;
  ib_euler_beam_bcs_t bcs;

  bcs.end[0] = IB_EULER_BEAM_BC_LEFT;
  bcs.end[1] = IB_EULER_BEAM_BC_RIGHT;
  bcs.type[0] = IB_EULER_BEAM_BC_FREE;
  bcs.type[1] = IB_EULER_BEAM_BC_FREE;

  requested_margin =
    fibre_margin > 0. ? fibre_margin : 0.5 * fibre_length + 2. * h_ibm;
  margin = min(requested_margin, 0.45 * L_fluid);

  rng_state = (unsigned)(fibre_seed ? fibre_seed : 1);

  for (int n = 0; n < fibre_count; n++) {
    coord center = {
      uniform_range(-0.5 * L_fluid + margin,
                     0.5 * L_fluid - margin),
      uniform_range(-0.5 * L_fluid + margin,
                     0.5 * L_fluid - margin),
      uniform_range(-0.5 * L_fluid + margin,
                     0.5 * L_fluid - margin)
    };
    const double azimuth = uniform_range(0., 2. * pi);
    const double direction_z = uniform_range(-1., 1.);
    const double direction_xy = sqrt(max(0., 1. - sq(direction_z)));
    const coord direction = {
      direction_xy * cos(azimuth),
      direction_xy * sin(azimuth),
      direction_z
    };
    const coord velocity = taylor_velocity(center);
    const coord directional_gradient =
      taylor_directional_velocity_gradient(center, direction);
    const double tangent_gradient =
      direction.x * directional_gradient.x +
      direction.y * directional_gradient.y +
      direction.z * directional_gradient.z;
    const coord compatible_gradient = {
      directional_gradient.x - tangent_gradient * direction.x,
      directional_gradient.y - tangent_gradient * direction.y,
      directional_gradient.z - tangent_gradient * direction.z
    };
    const coord start = {
      center.x - 0.5 * fibre_length * direction.x,
      center.y - 0.5 * fibre_length * direction.y,
      center.z - 0.5 * fibre_length * direction.z
    };

    int mesh_id = ibmeshmanager_add_mesh();
    IBMeshModel fibre_model =
      elff_euler_beam_huang_new_direction(
        fibre_length, fibre_bending_stiffness, fibre_linear_density,
        fibre_nodes, direction, bcs, start, fibre_pid);
    ib_euler_beam_huang_set_implicit_bending(
      (ib_euler_beam_huang_t)fibre_model.ctx,
      experiment_huang_implicit_bending);

    vertex_t initial_velocity[fibre_nodes];
    for (int ni = 0; ni < fibre_nodes; ++ni) {
      const double s = fibre_length * ni / (fibre_nodes - 1);
      const double offset = 0.5 * fibre_length - s;
      initial_velocity[ni] = (vertex_t){
        velocity.x + offset * compatible_gradient.x,
        velocity.y + offset * compatible_gradient.y,
        velocity.z + offset * compatible_gradient.z
      };
    }
    ib_euler_beam_huang_set_initial_velocity(
      (ib_euler_beam_huang_t)fibre_model.ctx,
      initial_velocity, fibre_nodes);
    ibmeshmanager_set_model(mesh_id, fibre_model);
  }

  foreach_ibnode_per_ibmesh() {
    mesh->depth = ibmlevel;
    node->depth = ibmlevel;
    foreach_dimension()
      ibval(gravity.x) = 0.;
  }
}

int main(int argc, char **argv) {
  input_file_register_option("basilisk.fluid", U0, PARAM_VALUE_DOUBLE);
  input_file_register_option("basilisk.fluid", Reynolds, PARAM_VALUE_DOUBLE);
  input_file_register_option("basilisk.fluid", dt_fluid, PARAM_VALUE_DOUBLE);
  input_file_register_option("basilisk.fluid", t_end, PARAM_VALUE_DOUBLE);
  input_file_register_option("basilisk.fluid", minlevel, PARAM_VALUE_INT);
  input_file_register_option("basilisk.fluid", maxlevel, PARAM_VALUE_INT);
  input_file_register_option("basilisk.fluid", ibmlevel, PARAM_VALUE_INT);
  input_file_register_option("basilisk.output", base_path, PARAM_VALUE_STRING);
  input_file_register_option("basilisk.output", dt_output, PARAM_VALUE_DOUBLE);
  input_file_register_option("fibre", fibre_count, PARAM_VALUE_INT);
  input_file_register_option("fibre", fibre_seed, PARAM_VALUE_INT);
  input_file_register_option("fibre", fibre_length, PARAM_VALUE_DOUBLE);
  input_file_register_option("fibre", fibre_aspect_ratio,
                             PARAM_VALUE_DOUBLE);
  input_file_register_option("fibre", fibre_nL3, PARAM_VALUE_DOUBLE);
  input_file_register_option("fibre", fibre_marker_spacing_ratio,
                             PARAM_VALUE_DOUBLE);
  input_file_register_option("fibre", fibre_density_ratio,
                             PARAM_VALUE_DOUBLE);
  input_file_register_option("fibre",
                             fibre_dimensionless_bending_rigidity,
                             PARAM_VALUE_DOUBLE);
  input_file_register_option("fibre", fibre_margin, PARAM_VALUE_DOUBLE);
  input_file_register_option("experiment", experiment_ib_force_relaxation,
                             PARAM_VALUE_DOUBLE);
  input_file_register_option("experiment", experiment_ib_richardson_iters,
                             PARAM_VALUE_INT);
  input_file_register_option("experiment", experiment_huang_implicit_bending,
                             PARAM_VALUE_INT);
  input_file_register_option("experiment",
                             experiment_adapt_velocity_tolerance,
                             PARAM_VALUE_DOUBLE);
  input_file_register_option("experiment", stats_interval, PARAM_VALUE_INT);

  int input_file_parse_result = input_file_parse_cli(argc, argv);
  if (input_file_parse_result != 0)
    return input_file_parse_result > 0 ? 0 : input_file_parse_result;

  input_file_print_options();
  input_file_apply_options();

  if (U0 <= 0. || Reynolds <= 0. || dt_fluid <= 0. ||
      t_end <= 0. || dt_output <= 0.) {
    fprintf(stderr, "fluid scales and time controls must be positive\n");
    return 2;
  }
  if (fibre_count <= 0 || fibre_nL3 <= 0. || fibre_length <= 0. ||
      fibre_aspect_ratio <= 0. || fibre_marker_spacing_ratio <= 0. ||
      fibre_dimensionless_bending_rigidity <= 0.) {
    fprintf(stderr, "fibre geometry and spacing must be positive\n");
    return 2;
  }
  if (maxlevel < ibmlevel) {
    fprintf(stderr, "maxlevel must be at least ibmlevel\n");
    return 2;
  }
  if (fibre_nodes <= 5) {
    fprintf(stderr,
            "fibre spacing yields only %d nodes; increase ibmlevel or "
            "fibre_length\n",
            fibre_nodes);
    return 2;
  }
  if (fibre_density_ratio <= 1.) {
    fprintf(stderr, "fibre_density_ratio must be greater than one\n");
    return 2;
  }
  if (experiment_adapt_velocity_tolerance <= 0.) {
    fprintf(stderr,
            "experiment_adapt_velocity_tolerance must be positive\n");
    return 2;
  }

  L0 = L_fluid * 1. [1];
  origin(-0.5 * L0, -0.5 * L0, -0.5 * L0);
  foreach_dimension()
    periodic(right);

  N = 1 << minlevel;
  DT = simulation_dt;
  mu = muv;
  ib_force_relaxation = experiment_ib_force_relaxation;
  ib_richardson_iters = experiment_ib_richardson_iters;
  display_control(Reynolds, 10, 1000);

  if (pid() == 0)
    fprintf(stderr,
            "nondimensional setup: box/Lf=%g L_TG/Lf=%g Re_TG=%g "
            "density_ratio=%g q=%g K_B=%g mu=%g EI=%g nL3=%g "
            "volume_fraction=%g count=%d nodes=%d ds/dx=%g "
            "nu=%g dt=%g t_end=%g dt_output=%g "
            "dt/T_TG=%g t_end/T_TG=%g dt_output/T_TG=%g\n",
            L_fluid / fibre_length, taylor_length_scale / fibre_length,
            Reynolds, fibre_density_ratio,
            fibre_line_density_ratio,
            fibre_dimensionless_bending_rigidity,
            fibre_linear_density, fibre_bending_stiffness,
            fibre_nL3, fibre_volume_fraction, fibre_count,
            fibre_nodes, fibre_ds / h_ibm, fluid_nu, simulation_dt,
            simulation_end_time, simulation_output_interval,
            dt_fluid, t_end, dt_output);

  run();
}

event properties(i++) {
  foreach_face()
    muv.x[] = fm.x[] * fluid_dynamic_viscosity;
}

event init(i = 0) {
  add_random_fibres();

  if (!restore_handler(base_path)) {
    set_taylor_green_field();
#if TREE
    adapt_wavelet_ibm(NULL, NULL, 0, 1, all, true);
#endif
  }
}

event logfile(i += stats_interval; t <= simulation_end_time) {
  double kinetic_energy = 0.;
  double enstrophy = 0.;
  double ib_force_x = 0.;
  double ib_force_y = 0.;
  double ib_force_z = 0.;

  foreach(reduction(+:kinetic_energy) reduction(+:enstrophy)
          reduction(+:ib_force_x) reduction(+:ib_force_y)
          reduction(+:ib_force_z)) {
    const double omega_x =
      (u.z[0,1] - u.z[0,-1] - u.y[0,0,1] + u.y[0,0,-1]) /
      (2. * Delta);
    const double omega_y =
      (u.x[0,0,1] - u.x[0,0,-1] - u.z[1] + u.z[-1]) /
      (2. * Delta);
    const double omega_z =
      (u.y[1] - u.y[-1] - u.x[0,1] + u.x[0,-1]) /
      (2. * Delta);
    kinetic_energy +=
      0.5 * (sq(u.x[]) + sq(u.y[]) + sq(u.z[])) * dv();
    enstrophy +=
      0.5 * (sq(omega_x) + sq(omega_y) + sq(omega_z)) * dv();
    ib_force_x += ibmf.x[] * dv();
    ib_force_y += ibmf.y[] * dv();
    ib_force_z += ibmf.z[] * dv();
  }

  if (pid() == 0)
    fprintf(stderr,
            "diag i=%d t/T_TG=%g cells=%ld fibres=%d nodes=%zu ke=%g "
            "enstrophy=%g ib_force=(%g,%g,%g)\n",
            i, t / taylor_time_scale, grid->tn, fibre_count,
            ibmm.pool.active.size, kinetic_energy,
            enstrophy, ib_force_x, ib_force_y, ib_force_z);
}

event output(i += 1; t <= simulation_end_time) {
  scalar l2[];
  vector omega[];

  foreach() {
    omega.x[] =
      (u.z[0,1] - u.z[0,-1] - u.y[0,0,1] + u.y[0,0,-1]) /
      (2. * Delta);
    omega.y[] =
      (u.x[0,0,1] - u.x[0,0,-1] - u.z[1] + u.z[-1]) /
      (2. * Delta);
    omega.z[] =
      (u.y[1] - u.y[-1] - u.x[0,1] + u.x[0,-1]) /
      (2. * Delta);
  }
  boundary((scalar *){omega});
  lambda2(u, l2);

#if TREE
  output_hdf_htg({l2, p}, {u, omega, ibmf}, base_path);
#else
  output_hdf_imagedata({l2, p}, {u, omega, ibmf}, base_path);
#endif
  output_hdf_pd(NULL, (IBvector[]){eulvel, nforce, nvel, {{-1}}}, base_path);
}

#if TREE
event adapt(i++) {
  const double tol = experiment_adapt_velocity_tolerance * U0;
  adapt_wavelet_ibm({u}, (double[]){tol, tol, tol}, maxlevel, minlevel);
}
#endif

event checkpoint_event(i++, last) {
  return checkpoint_handler(t, i, base_path);
}

event stop(t = simulation_end_time) {
  if (pid() == 0)
    fprintf(stderr, "finished decaying Taylor vortex fibres at t=%g i=%d\n",
            t / taylor_time_scale, i);
  return 1;
}
