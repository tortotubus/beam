/* Two-dimensional analogue of src/examples/isotropic.c: periodic flow with
 * mean-subtracted linear forcing and a free Huang fibre. 2D dynamics and a
 * single oriented fibre do not reproduce the 3D isotropic benchmark. */
#include "library/ibm/IBMeshManager.h"
#include "library/ibm/navier-stokes/centered-mdf.h"
#include "library/elff/elff.h"
#include "library/ibm/IBOutput.h"
#include "library/io/output-vtk.h"
#include "library/io/params/params-cli.h"

#if dimension != 2 || !TREE
#error isotropic_fibre requires a two-dimensional tree grid
#endif

double L_fluid = 2.*pi;
double dt_fluid = 0.001;
double fluid_nu = 0.01;
double U0 = 1.;
double forcing_rate = 0.1;
double velocity_tolerance = 0.01;
int minlevel = 5;
int maxlevel = 8;
int ibmlevel = 8;
double huang_gamma = 0.001;
double huang_rho = 1.5;
double b_length = 1.;
/* Match flapping_fibre: (1/64)/(16/2048) = 2.
 * Zero selects the nearest integer node count for the requested spacing. */
int b_nodes = 0;
double ds_dx = 2.;
double b_theta = pi/4.;
double t_end = 10.;
double output_interval = 0.1;
int stats_interval = 10;
double force_relaxation = 0.1;
#if _MPI
char * base_path = "isotropic_fibre_mpi_output";
#else
char * base_path = "isotropic_fibre_serial_output";
#endif

face vector muv[], av[];

int main (int argc, char ** argv) {
  input_file_register_option("basilisk.fluid", L_fluid, PARAM_VALUE_DOUBLE);
  input_file_register_option("basilisk.fluid", dt_fluid, PARAM_VALUE_DOUBLE);
  input_file_register_option("basilisk.fluid", fluid_nu, PARAM_VALUE_DOUBLE);
  input_file_register_option("basilisk.fluid", U0, PARAM_VALUE_DOUBLE);
  input_file_register_option("basilisk.fluid", forcing_rate, PARAM_VALUE_DOUBLE);
  input_file_register_option("basilisk.fluid", velocity_tolerance, PARAM_VALUE_DOUBLE);
  input_file_register_option("basilisk.fluid", minlevel, PARAM_VALUE_INT);
  input_file_register_option("basilisk.fluid", maxlevel, PARAM_VALUE_INT);
  input_file_register_option("basilisk.fluid", ibmlevel, PARAM_VALUE_INT);
  input_file_register_option("huang", huang_gamma, PARAM_VALUE_DOUBLE);
  input_file_register_option("huang", huang_rho, PARAM_VALUE_DOUBLE);
  input_file_register_option("beam", b_length, PARAM_VALUE_DOUBLE);
  input_file_register_option("beam", b_nodes, PARAM_VALUE_INT);
  input_file_register_option("beam", ds_dx, PARAM_VALUE_DOUBLE);
  input_file_register_option("beam", b_theta, PARAM_VALUE_DOUBLE);
  input_file_register_option("basilisk.output", base_path, PARAM_VALUE_STRING);
  input_file_register_option("basilisk.output", output_interval, PARAM_VALUE_DOUBLE);
  input_file_register_option("basilisk.output", stats_interval, PARAM_VALUE_INT);
  input_file_register_option("experiment", t_end, PARAM_VALUE_DOUBLE);
  input_file_register_option("beam", force_relaxation, PARAM_VALUE_DOUBLE);
  int result = input_file_parse_cli(argc, argv);
  if (result != 0)
    return result > 0 ? 0 : result;
  input_file_apply_options();
  input_file_print_options();

  if (!(L_fluid > 0. && dt_fluid > 0. && fluid_nu > 0. && U0 >= 0. &&
        forcing_rate >= 0. && velocity_tolerance > 0. &&
        huang_gamma >= 0. && huang_rho > 0. && b_length > 0. &&
        b_length < L_fluid && force_relaxation > 0. && force_relaxation <= 1. && (b_nodes == 0 || b_nodes >= 5) && ds_dx > 0. && isfinite(ds_dx) && isfinite(b_theta) &&
        minlevel >= 2 && maxlevel >= minlevel && maxlevel <= 12 &&
        ibmlevel >= minlevel && ibmlevel <= 12 && t_end > 0. &&
        output_interval > 0. && stats_interval > 0. && base_path && base_path[0])) {
    if (pid() == 0)
      fprintf(stderr, "Invalid isotropic_fibre parameters (tree levels must be 2..12).\n");
    return 1;
  }
  if (b_nodes == 0) {
    double segments = b_length/(ds_dx*L_fluid/(1 << ibmlevel));
    if (!isfinite(segments) || segments < 3.5 || segments > 1000000.) {
      if (pid() == 0)
        fprintf(stderr, "Requested ds/dx needs 4..1000000 segments; adjust ibmlevel or ds_dx.\n");
      return 1;
    }
    b_nodes = (int) floor(segments + 0.5) + 1;
  }
  if (pid() == 0)
    fprintf(stderr, "Fibre nodes=%d ds/dx=%g (target=%g)\n", b_nodes,
            b_length/(b_nodes - 1)/(L_fluid/(1 << ibmlevel)), ds_dx);
  L0 = L_fluid;
  origin(-L0/2., -L0/2.);
  foreach_dimension()
    periodic(right);
  N = 1 << minlevel;
  DT = dt_fluid;
  ib_force_relaxation = force_relaxation;
  mu = muv;
  a = av;
  run();
}

event properties (i++) {
  foreach_face()
    muv.x[] = fm.x[]*fluid_nu;
}

/* Analytic curl of a periodic streamfunction. Equal amplitudes for rotated
 * wavevectors avoid a preferred coordinate axis; multiple shells allow nonlinear
 * interactions, unlike a single Taylor--Green mode. Phases are deterministic
 * so serial and MPI runs start from the same field. U0 is RMS speed. */
void initial_velocity () {
  int kx[8] = {1, 0, 1, 1, 2, 1, 2, 1};
  int ky[8] = {0, 1, 1, -1, 1, 2, -1, -2};
  double phase[8] = {0.2, 1.1, 2.4, 0.7, 1.7, 2.8, 0.4, 2.1};
  foreach() {
    u.x[] = u.y[] = 0.;
    for (int m = 0; m < 8; m++) {
      double q = 2.*pi*(kx[m]*x + ky[m]*y)/L0 + phase[m];
      double amplitude = U0/2./sqrt(sq(kx[m]) + sq(ky[m]));
      u.x[] += amplitude*ky[m]*cos(q);
      u.y[] -= amplitude*kx[m]*cos(q);
    }
  }
  boundary((scalar *){u});
}

event init (i = 0) {
  int m_id = ibmeshmanager_add_mesh();
  ib_euler_beam_bcs_t bcs = elff_euler_beam_bcs_types(
    IB_EULER_BEAM_BC_FREE, IB_EULER_BEAM_BC_FREE);
  coord start = {-0.5*b_length*cos(b_theta), -0.5*b_length*sin(b_theta)};
  IBMeshModel model = elff_euler_beam_huang_new_theta(
    b_length, huang_gamma, huang_rho, b_nodes, b_theta, bcs, start, 0);
  ib_euler_beam_huang_set_implicit_bending((ib_euler_beam_huang_t)model.ctx, 1);
  ibmeshmanager_set_model(m_id, model);
  foreach_ibnode_per_ibmesh() {
    mesh->depth = node->depth = ibmlevel;
    foreach_dimension()
      ibval(gravity.x) = 0.;
  }
  if (!restore_handler(base_path)) {
    initial_velocity();
    adapt_wavelet_ibm({u}, (double[]){velocity_tolerance, velocity_tolerance},
                      maxlevel, minlevel, all, true);
    initial_velocity();
  }
}

event acceleration (i++) {
  coord ubar;
  foreach_dimension() {
    stats s = statsf(u.x);
    ubar.x = s.sum/s.volume;
  }
  foreach_face()
    av.x[] += forcing_rate*((u.x[] + u.x[-1])/2. - ubar.x);
}

event statistics (i += stats_interval) {
  coord ubar;
  foreach_dimension() {
    stats s = statsf(u.x);
    ubar.x = s.sum/s.volume;
  }
  scalar omega[];
  vorticity(u, omega);
  double ke = 0., dissipation = 0., enstrophy = 0., volume = 0.;
  foreach(reduction(+:ke) reduction(+:dissipation)
          reduction(+:enstrophy) reduction(+:volume)) {
    volume += dv();
    enstrophy += dv()*sq(omega[]);
    foreach_dimension() {
      ke += dv()*sq(u.x[] - ubar.x);
      dissipation += dv()*(sq(u.x[1] - u.x[-1]) +
                          sq(u.x[0,1] - u.x[0,-1]))/sq(2.*Delta);
    }
  }
  ke /= 2.*volume;
  dissipation *= fluid_nu/volume;
  enstrophy /= 2.*volume;
  if (pid() == 0) {
    create_path(base_path);
    char fname[4096];
    snprintf(fname, sizeof(fname), "%s/flow.csv", base_path);
    FILE * fp = fopen(fname, "a+");
    if (!fp) {
      perror(fname);
      return 1;
    }
    fseek(fp, 0, SEEK_END);
    if (ftell(fp) == 0)
      fprintf(fp, "i,t,kinetic_energy,dissipation,enstrophy,forcing_power,ux_mean,uy_mean\n");
    fprintf(fp, "%d,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g\n",
            i, t, ke, dissipation, enstrophy, 2.*forcing_rate*ke, ubar.x, ubar.y);
    fclose(fp);
    fprintf(stderr, "i=%d t=%g energy=%g dissipation=%g\n", i, t, ke, dissipation);
  }
}

event output (t += output_interval) {
  scalar omega[];
  vorticity(u, omega);
  output_hdf_htg({omega, p}, {u, ibmf}, base_path);
  output_hdf_pd(NULL, (IBvector[]){eulvel, nforce, nvel, {{-1}}}, base_path);
}

event adapt (i++) {
  adapt_wavelet_ibm({u}, (double[]){velocity_tolerance, velocity_tolerance},
                    maxlevel, minlevel);
}

event checkpoint_event (i++, last) {
  return checkpoint_handler(t, i, base_path);
}

event stop (t = t_end) {
  return 1;
}
