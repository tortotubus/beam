/**
 * Fluid-only Eswaran-Pope forcing in a triply periodic cube.
 * The historical filename is retained; no ELFF or immersed-boundary code is used.
 * Usage: executable [minlevel [end_time [sigma [correlation_time [seed [maxlevel [output_interval]]]]]]]
 * Output: volume-weighted diagnostics in energy.csv (overwritten on a fresh run).
 */
#include "navier-stokes/centered.h"
#include "isotropic_forcing.h"
#include "library/io/output-vtk-htg.h"
#include "lambda2.h"
#include <errno.h>
#include <limits.h>

#if dimension != 3
#error "Eswaran-Pope forcing requires a three-dimensional grid"
#endif

static int minlevel = 9, maxlevel = 11;
static double end_time = 10., output_interval = 0.05;
static double forcing_sigma = 0.05, forcing_time = 1.;
static uint64_t forcing_seed = 1;
static double viscosity_value = 0.01, adapt_tolerance = 0.01;
static EPForcing forcing;
static FILE *energy_file;
face vector acceleration_field[], stochastic_acceleration[];
vector velocity_before_projection[], forcing_at_cells[];

/** @brief Broadcast coefficients, never grid-dependent random draws. */
static void broadcast_forcing (void) {
#if _MPI
  MPI_Bcast(forcing.b, 54, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&forcing.counter, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
#endif
  ep_project(&forcing);
}

/** @brief Parse finite numeric arguments without silent truncation. */
static double numeric_argument (const char *text) {
  char *end;
  errno = 0;
  double value = strtod(text, &end);
  if (errno || end == text || *end || !isfinite(value)) {
    fprintf(stderr, "Invalid numeric argument: %s\n", text);
    exit(2);
  }
  return value;
}

int main (int argc, char **argv) {
  if (argc > 8) {
    if (pid() == 0) fprintf(stderr,
      "Usage: %s [minlevel [end_time [sigma [correlation_time [seed [maxlevel [output_interval]]]]]]]\n", argv[0]);
    return 2;
  }
  double level = argc > 1 ? numeric_argument(argv[1]) : minlevel;
  double maximum = argc > 6 ? numeric_argument(argv[6]) : level;
  if (level < 3 || level > 12 || floor(level) != level ||
      maximum < level || maximum > 12 || floor(maximum) != maximum) {
    if (pid() == 0) fprintf(stderr, "Require integer 3 <= minlevel <= maxlevel <= 12.\n");
    return 2;
  }
  minlevel = (int)level; maxlevel = (int)maximum;
  if (argc > 2) end_time = numeric_argument(argv[2]);
  if (argc > 3) forcing_sigma = numeric_argument(argv[3]);
  if (argc > 4) forcing_time = numeric_argument(argv[4]);
  if (argc > 7) output_interval = numeric_argument(argv[7]);
  if (argc > 5) {
    char *end;
    errno = 0;
    unsigned long long seed = strtoull(argv[5], &end, 10);
    if (errno || end == argv[5] || *end || argv[5][0] == '-') {
      if (pid() == 0) fprintf(stderr, "Seed must be an unsigned integer.\n");
      return 2;
    }
    forcing_seed = (uint64_t)seed;
  }
  L0 = 2.*pi;
  origin(-L0/2., -L0/2., -L0/2.);
  foreach_dimension() periodic(right);
  N = 1 << minlevel;
  DT = 0.001;
  mu[] = {viscosity_value, viscosity_value, viscosity_value};
  a = acceleration_field;
  if (end_time <= 0. || output_interval <= 0. || ep_init(&forcing, L0, forcing_time, forcing_sigma,
                               forcing_seed, pid() == 0)) {
    if (pid() == 0) fprintf(stderr, "Require end_time > 0, output_interval > 0, correlation_time > 0 and sigma >= 0.\n");
    return 2;
  }
  broadcast_forcing();
  run();
  return 0;
}

/** @brief Start the fluid at rest; the OU coefficients start stationary. */
event init (i = 0) {
  if (pid() == 0) {
    energy_file = fopen("energy.csv", "w");
    if (!energy_file) perror("energy.csv");
  }
  int opened = pid() == 0 ? energy_file != NULL : 0;
#if _MPI
  MPI_Bcast(&opened, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
  if (!opened) exit(1);
  if (pid() == 0) {
    fprintf(energy_file, "i,t,dt,K,P_turb_midpoint,epsilon,force_rms,force_div_rms\n");
    fprintf(stderr, "Eswaran-Pope: sigma=%g Tf=%g seed=%llu expected_force_rms=%g\n",
      forcing_sigma, forcing_time, (unsigned long long)forcing_seed, 6.*sqrt(2.)*forcing_sigma);
  }
  foreach() foreach_dimension() u.x[] = 0.;
}

/** @brief Update once per fluid step, before centered.h consumes acceleration.
 * Each explicit face loop uses Basilisk's physical face-centre coordinates.
 * centered.h clears the nonconstant a field before this inherited event.
 */
event acceleration (i++) {
  if (pid() == 0 && ep_advance(&forcing, dt)) {
    fprintf(stderr, "Invalid OU timestep\n"); exit(1);
  }
  broadcast_forcing();
  foreach_face(x)
    stochastic_acceleration.x[] = ep_component(&forcing, x-X0, y-Y0, z-Z0, 0);
  foreach_face(y)
    stochastic_acceleration.y[] = ep_component(&forcing, x-X0, y-Y0, z-Z0, 1);
  foreach_face(z)
    stochastic_acceleration.z[] = ep_component(&forcing, x-X0, y-Y0, z-Z0, 2);
  boundary((scalar *){stochastic_acceleration});
  foreach_face() a.x[] += stochastic_acceleration.x[];
  foreach() foreach_dimension() {
    velocity_before_projection.x[] = u.x[];
    forcing_at_cells.x[] = (fm.x[]*stochastic_acceleration.x[] +
      fm.x[1]*stochastic_acceleration.x[1])/(fm.x[] + fm.x[1] + SEPS);
  }
}

/** @brief Measure the resolved work, strain dissipation and discrete divergence.
 * Midpoint velocity straddles the final pressure/acceleration correction.
 * Work uses the same metric-weighted face-to-cell reconstruction as g; this
 * is a diagnostic estimate, not an exact discrete energy-budget identity.
 */
event end_timestep (i++) {
  boundary((scalar *){u});
  double energy = 0., dissipation = 0., power = 0., volume = 0.;
  double force2 = 0., div2 = 0.;
  foreach(reduction(+:energy) reduction(+:dissipation) reduction(+:power)
          reduction(+:volume) reduction(+:force2) reduction(+:div2)) {
    double vol = dv(), divergence = 0.;
    volume += vol;
    foreach_dimension() {
      energy += vol*sq(u.x[])/2.;
      power += vol*0.5*(u.x[] + velocity_before_projection.x[])*forcing_at_cells.x[];
      force2 += vol*sq(forcing_at_cells.x[]);
      divergence += (stochastic_acceleration.x[1] - stochastic_acceleration.x[])/Delta;
    }
    div2 += vol*sq(divergence);
    double sxx = (u.x[1] - u.x[-1])/(2.*Delta);
    double syy = (u.y[0,1] - u.y[0,-1])/(2.*Delta);
    double szz = (u.z[0,0,1] - u.z[0,0,-1])/(2.*Delta);
    double sxy = (u.x[0,1]-u.x[0,-1] + u.y[1]-u.y[-1])/(4.*Delta);
    double sxz = (u.x[0,0,1]-u.x[0,0,-1] + u.z[1]-u.z[-1])/(4.*Delta);
    double syz = (u.y[0,0,1]-u.y[0,0,-1] + u.z[0,1]-u.z[0,-1])/(4.*Delta);
    dissipation += vol*2.*viscosity_value*(sq(sxx)+sq(syy)+sq(szz)+2.*(sq(sxy)+sq(sxz)+sq(syz)));
  }
  if (pid() == 0) {
    fprintf(energy_file, "%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
      i, t+dt, dt, energy/volume, power/volume, dissipation/volume,
      sqrt(force2/volume), sqrt(div2/volume));
    fflush(energy_file);
    if (i % 100 == 0)
      fprintf(stderr, "HIT i=%d t=%g K=%g P=%g epsilon=%g\n",
        i, t+dt, energy/volume, power/volume, dissipation/volume);
  }
}

#if TREE
event adapt (i++) {
  if (maxlevel > minlevel)
    adapt_wavelet((scalar *){u}, (double[]){adapt_tolerance, adapt_tolerance, adapt_tolerance},
                  maxlevel, minlevel);
}
#endif

/** @brief ParaView HTG time series of pressure, velocity and vortex structure. */
event output (t = 0; t += output_interval; t <= end_time) {
  scalar l2[];
  lambda2(u, l2);
  output_hdf_htg({p, l2}, {u}, "isotropic_flow_output");
}

event stop (t = end_time) { return 1; }
event cleanup_output (t = end) {
  if (energy_file) { fclose(energy_file); energy_file = NULL; }
}
