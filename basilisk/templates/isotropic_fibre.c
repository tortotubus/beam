#include "library/ibm/IBMeshManager.h"
#include "library/ibm/navier-stokes/centered-mdf.h"
#include "library/elff/elff.h"
#include "library/ibm/IBOutput.h"
#include "library/io/output-vtk.h"
#include "library/io/params/params-cli.h"

#include "lambda2.h"
 

double dt_fluid = 0.001;
double L_fluid = 2. * pi;
double fluid_nu = 1.59e-3;

int minlevel = 8;
int maxlevel = 8;
int ibmlevel = 8;

/* Direct beam parameters, retained while fibre creation is disabled. */
double b_length = 1.0;
int b_nodes = 21;
double b_theta = 0.0 * pi;
double b_EI = 1e-4;
double b_mu = 1 + 1e-3;
double b_gravity = 0.;

double turbulence_u0 = 0.514;
double turbulence_forcing_rate = 0.056776;
double turbulence_adapt_tolerance = 5e-3;

double experiment_output_interval = 0.05;
double experiment_t_end = 50.;

char *base_path = "isotropic_flow_output";

#define b_ds (b_length / (b_nodes - 1))
#define turbulence_k (2. * pi / L_fluid)

face vector muv[];
face vector turbulence_acceleration[];

static coord isotropic_base_velocity(coord position) {
  const double k = turbulence_k;
  return (coord){turbulence_u0 * (cos(k * position.y) + sin(k * position.z)),
                 turbulence_u0 * (sin(k * position.x) + cos(k * position.z)),
                 turbulence_u0 * (cos(k * position.x) + sin(k * position.y))};
}

static void add_isotropic_fibre(void) {
  const coord center = {0.};
  const coord direction = {cos(b_theta), sin(b_theta), 0.};
  const coord start = {-0.5 * b_length * direction.x,
                       -0.5 * b_length * direction.y,
                       -0.5 * b_length * direction.z};

  ib_euler_beam_bcs_t bcs = {0};
  bcs.end[0] = IB_EULER_BEAM_BC_LEFT;
  bcs.end[1] = IB_EULER_BEAM_BC_RIGHT;
  bcs.type[0] = IB_EULER_BEAM_BC_FREE;
  bcs.type[1] = IB_EULER_BEAM_BC_FREE;

  IBMeshModel model = elff_euler_beam_huang_new_direction(
      b_length, b_EI, b_mu, b_nodes, direction, bcs, start, 0);

  /* Project out tangential strain from the initial nodal velocities. */
  const coord center_velocity = isotropic_base_velocity(center);
  const double center_tangent = center_velocity.x * direction.x +
                                center_velocity.y * direction.y +
                                center_velocity.z * direction.z;
  vertex_t initial_velocity[b_nodes];

  for (int ni = 0; ni < b_nodes; ++ni) {
    const double s = ni * b_ds;
    const coord position = {start.x + s * direction.x,
                            start.y + s * direction.y,
                            start.z + s * direction.z};
    coord velocity = isotropic_base_velocity(position);
    const double tangent = velocity.x * direction.x + velocity.y * direction.y +
                           velocity.z * direction.z;
    const double correction = center_tangent - tangent;
    velocity.x += correction * direction.x;
    velocity.y += correction * direction.y;
    velocity.z += correction * direction.z;
    initial_velocity[ni] = (vertex_t){velocity.x, velocity.y, velocity.z};
  }

  ib_euler_beam_huang_set_initial_velocity((ib_euler_beam_huang_t)model.ctx,
                                           initial_velocity, b_nodes);

  const int mesh_id = ibmeshmanager_add_mesh();
  ibmeshmanager_set_model(mesh_id, model);

  foreach_ibnode_per_ibmesh() {
    mesh->depth = ibmlevel;
    node->depth = ibmlevel;
    foreach_dimension() ibval(gravity.x) = 0.;
    ibval(gravity.y) = b_gravity;
  }
}

int main(int argc, char **argv) {
  input_file_register_option("basilisk.fluid", L_fluid, PARAM_VALUE_DOUBLE);
  input_file_register_option("basilisk.fluid", dt_fluid, PARAM_VALUE_DOUBLE);
  input_file_register_option("basilisk.fluid", fluid_nu, PARAM_VALUE_DOUBLE);
  input_file_register_option("basilisk.fluid", minlevel, PARAM_VALUE_INT);
  input_file_register_option("basilisk.fluid", maxlevel, PARAM_VALUE_INT);
  input_file_register_option("basilisk.fluid", ibmlevel, PARAM_VALUE_INT);
  input_file_register_option("basilisk.output", base_path, PARAM_VALUE_STRING);

  input_file_register_option("beam", b_length, PARAM_VALUE_DOUBLE);
  input_file_register_option("beam", b_nodes, PARAM_VALUE_INT);
  input_file_register_option("beam", b_theta, PARAM_VALUE_DOUBLE);
  input_file_register_option("beam", b_EI, PARAM_VALUE_DOUBLE);
  input_file_register_option("beam", b_mu, PARAM_VALUE_DOUBLE);
  input_file_register_option("beam", b_gravity, PARAM_VALUE_DOUBLE);

  input_file_register_option("turbulence", turbulence_u0, PARAM_VALUE_DOUBLE);
  input_file_register_option("turbulence", turbulence_forcing_rate,
                             PARAM_VALUE_DOUBLE);
  input_file_register_option("turbulence", turbulence_adapt_tolerance,
                             PARAM_VALUE_DOUBLE);

  input_file_register_option("experiment", experiment_output_interval,
                             PARAM_VALUE_DOUBLE);
  input_file_register_option("experiment", experiment_t_end,
                             PARAM_VALUE_DOUBLE);

  const int parse_result = input_file_parse_cli(argc, argv);
  if (parse_result != 0)
    return parse_result > 0 ? 0 : parse_result;

  input_file_print_options();
  input_file_apply_options();

  if (L_fluid <= 0. || dt_fluid <= 0. || fluid_nu <= 0. || b_length <= 0. ||
      b_nodes <= 5 || b_EI <= 0. || b_mu <= 0. || turbulence_u0 <= 0. ||
      turbulence_forcing_rate < 0. || turbulence_adapt_tolerance <= 0. ||
      experiment_output_interval <= 0. || experiment_t_end <= 0.) {
    fprintf(stderr, "isotropic fibre parameters must be positive\n");
    return 2;
  }
  if (minlevel < 0 || maxlevel < minlevel || ibmlevel < minlevel) {
    fprintf(stderr, "levels must satisfy 0 <= minlevel <= maxlevel and "
                    "minlevel <= ibmlevel\n");
    return 2;
  }

  L0 = L_fluid;
#if dimension == 2
  origin(-L_fluid / 2., -L_fluid / 2.);
#else
  origin(-L_fluid / 2., -L_fluid / 2., -L_fluid / 2.);
#endif
  foreach_dimension() periodic(right);

  N = 1 << minlevel;
  DT = dt_fluid;
  mu = muv;
  a = turbulence_acceleration;

  if (pid() == 0)
    fprintf(stderr,
            "isotropic flow (fibre disabled): dimension=%d L0=%g u0=%g "
            "nu=%g forcing_rate=%g dt=%g levels=%d/%d t_end=%g\n"
            "retained beam: length=%g nodes=%d mu=%g EI=%g gravity=%g "
            "ibmlevel=%d\n",
            dimension, L_fluid, turbulence_u0, fluid_nu,
            turbulence_forcing_rate, dt_fluid, minlevel, maxlevel,
            experiment_t_end, b_length, b_nodes, b_mu, b_EI, b_gravity,
            ibmlevel);

  run();
}

event properties(i++) { foreach_face() muv.x[] = fm.x[] * fluid_nu; }

event init(i = 0) {
  add_isotropic_fibre();

  if (!restore_handler(base_path)) {
    foreach () {
      const coord velocity = isotropic_base_velocity((coord){x, y, z});
      foreach_dimension() u.x[] = velocity.x;
    }
    boundary((scalar *){u});
  }
}

event acceleration(i++) {
  coord mean_velocity = {0.};
  foreach_dimension() {
    const stats velocity_stats = statsf(u.x);
    mean_velocity.x = velocity_stats.sum / velocity_stats.volume;
  }

  foreach_face() turbulence_acceleration.x[] +=
      turbulence_forcing_rate * ((u.x[] + u.x[-1]) / 2. - mean_velocity.x);
}

event logfile(i++; t <= experiment_t_end) {
  coord mean_velocity = {0.};
  foreach_dimension() {
    const stats velocity_stats = statsf(u.x);
    mean_velocity.x = velocity_stats.sum / velocity_stats.volume;
  }

  double kinetic_energy = 0.;
  double dissipation = 0.;
  double volume = 0.;
  foreach (reduction(+ : kinetic_energy) reduction(+ : dissipation)
               reduction(+ : volume)) {
    const double cell_volume = dv();
    volume += cell_volume;
    foreach_dimension() {
      kinetic_energy += cell_volume * sq(u.x[] - mean_velocity.x);
      double gradient_squared = sq(u.x[1] - u.x[-1]);
      gradient_squared += sq(u.x[0, 1] - u.x[0, -1]);
#if dimension == 3
      gradient_squared += sq(u.x[0, 0, 1] - u.x[0, 0, -1]);
#endif
      dissipation += cell_volume * gradient_squared / sq(2. * Delta);
    }
  }
  kinetic_energy /= 2. * volume;
  dissipation *= fluid_nu / volume;

#if dimension == 3
  const double reynolds_lambda = dissipation > 0.
                                     ? 2. / 3. * kinetic_energy / fluid_nu *
                                           sqrt(15. * fluid_nu / dissipation)
                                     : 0.;
#else
  const double reynolds_lambda =
      dissipation > 0.
          ? kinetic_energy / fluid_nu * sqrt(4. * fluid_nu / dissipation)
          : 0.;
#endif

  if (pid() == 0) {
    if (i == 0)
      fprintf(stderr, "t dissipation energy Reynolds_lambda\n");
    fprintf(stderr, "%g %g %g %g\n", t, dissipation, kinetic_energy,
            reynolds_lambda);
  }
}

#if TREE
event adapt(i++) {
  if (maxlevel == minlevel)
    return 0;
#if dimension == 2
  adapt_wavelet(
      (scalar *){u},
      (double[]){turbulence_adapt_tolerance, turbulence_adapt_tolerance},
      maxlevel, minlevel);
#else
  adapt_wavelet((scalar *){u},
                (double[]){turbulence_adapt_tolerance,
                           turbulence_adapt_tolerance,
                           turbulence_adapt_tolerance},
                maxlevel, minlevel);
#endif
}
#endif

event output(i += 1; t <= experiment_t_end) {
  scalar omega_z[];
  vorticity(u, omega_z);
  scalar l2[];
  lambda2(u, l2);
#if TREE
  output_hdf_htg({l2, omega_z, p}, {u, ibmf}, base_path);
#else
  output_hdf_imagedata({l2, omega_z, p}, {u, ibmf}, base_path);
#endif
  output_hdf_pd(NULL, (IBvector[]){eulvel, nforce, nvel, {{-1}}}, base_path);
}

event checkpoint_event(i++, last) {
  return checkpoint_handler(t, i, base_path);
}

event stop(t = experiment_t_end) {
  if (pid() == 0)
    fprintf(stderr, "finished isotropic flow at t=%g i=%d\n", t, i);
  return 1;
}
