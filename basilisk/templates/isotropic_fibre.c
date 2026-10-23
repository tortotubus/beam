/**
 * Fluid-only Eswaran-Pope forcing in a triply periodic cube.
 * The historical filename is retained; no ELFF or immersed-boundary code is used.
 * Usage: executable [minlevel [end_time [sigma [correlation_time [seed [maxlevel [output_interval]]]]]]]
 * Output: volume-weighted diagnostics in energy.csv (overwritten on a fresh run).
 *
 * @par Validation and comparison references
 * Runtime checks below remain enabled with NDEBUG. They check algebra and
 * spatial sampling, not whether a flow is resolved or statistically stationary.
 * First run a fixed grid (minlevel == maxlevel); repeat at finer resolution
 * and smaller DT before comparing an adaptive run with that baseline.
 * A sigma=0 run from rest must retain zero velocity and forcing.
 *
 * The forcing reference is Eswaran & Pope, Computers & Fluids 16 (1988),
 * 257-278, https://doi.org/10.1016/0045-7930(88)90013-8 . Match the forced
 * wavenumber band, Fourier normalization, viscosity and correlation time
 * before comparing statistics; their plotted energies are not universal targets.
 * For this implementation each unprojected real OU component has
 * \f( E[b^2]=\sigma^2 \f) and correlation
 * \f( E[b(t)b(t+\tau)]=\sigma^2 e^{-|\tau|/T_f} \f).
 * See tests/unit/basilisk/EswaranPope.cpp for ensemble tests of these identities.
 * Nine conjugate pairs and a rank-two transverse projector give
 * \f( E[\langle |f|^2\rangle]=72\sigma^2 \f) in the continuum.
 * Compare a long-time/ensemble mean of force_rms squared, not the mean RMS,
 * with this value, allowing for face-to-cell filtering and correlated samples.
 * No instantaneous random realization is required to equal the ensemble mean.
 *
 * For the fluid, compare the periodic incompressible kinetic-energy balance
 * \f[ (K(t_1)-K(t_0))/(t_1-t_0)
 *       = \overline{P-\epsilon},\qquad
 *       \epsilon=2\nu\langle S_{ij}S_{ij}\rangle. \f]
 * Integrate the CSV rates using dt, discard the startup transient for stationary
 * averages, and check that the residual decreases under space/time refinement.
 * The reported work and strain are diagnostic approximations, so this is not
 * a roundoff-level assertion. At stationarity mean power balances dissipation.
 * Compare spectra and component statistics against an independently resolved
 * spectral DNS with matched forcing; passing these checks alone does not prove
 * isotropy or DNS resolution. Basilisk's src/examples/isotropic.c is a useful
 * solver comparison example, but uses linear forcing and different initial data.
 */
#include "navier-stokes/centered.h"
#include "isotropic_forcing.h"
#include "library/io/output-vtk-htg.h"
#include "lambda2.h"
#include <errno.h>
#include <limits.h>
#include <float.h>

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

/** @brief Fail on all MPI processes even when assertions are compiled out. */
static void validation_require (int valid, const char *message) {
  if (valid) return;
  fprintf(stderr, "Isotropic validation failed on rank %d: %s\n", pid(), message);
#if _MPI
  MPI_Abort(MPI_COMM_WORLD, 1);
#endif
  exit(1);
}

/** @brief Broadcast coefficients, never grid-dependent random draws. */
static void broadcast_forcing (void) {
#if _MPI
  MPI_Bcast(forcing.b, 54, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&forcing.counter, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
#endif
  ep_project(&forcing);
  for (int m = 0; m < 9; ++m)
    for (int part = 0; part < 2; ++part) {
      double dot = 0., scale = 0.;
      for (int d = 0; d < 3; ++d) {
        double b = forcing.b[m][3*part+d];
        double projected = forcing.transverse[m][3*part+d];
        validation_require(isfinite(b) && isfinite(projected),
                           "nonfinite Fourier coefficient");
        dot += ep_wave_numbers[m][d]*projected;
        scale += fabs(b);
      }
      validation_require(fabs(dot) <= 64.*DBL_EPSILON*scale,
                         "forcing is not transverse to its wavevector");
    }
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
  double force2 = 0., div2 = 0., mean_fx = 0., mean_fy = 0., mean_fz = 0.;
  foreach(reduction(+:energy) reduction(+:dissipation) reduction(+:power)
          reduction(+:volume) reduction(+:force2) reduction(+:div2)
          reduction(+:mean_fx) reduction(+:mean_fy) reduction(+:mean_fz)) {
    double vol = dv(), divergence = 0.;
    volume += vol;
    mean_fx += vol*forcing_at_cells.x[];
    mean_fy += vol*forcing_at_cells.y[];
    mean_fz += vol*forcing_at_cells.z[];
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
  validation_require(isfinite(volume) && volume > 0. &&
                     isfinite(energy) && energy >= 0. &&
                     isfinite(dissipation) && dissipation >= 0. &&
                     isfinite(power) && isfinite(force2) && force2 >= 0. &&
                     isfinite(div2) && div2 >= 0., "invalid flow diagnostics");
  if (forcing_sigma == 0.)
    validation_require(energy == 0. && force2 == 0. && dissipation == 0.,
                       "unforced fluid did not remain at rest");
  if (minlevel == maxlevel) {
    /** Parseval with the actual face-to-cell filter: averaging opposite faces
     * multiplies component d of each mode by cos(k_d h/2). Orthogonality is
     * exact on these uniform periodic grids (N >= 8); the factor two accounts
     * for the conjugate pair. This detects missing modes, amplitude errors and
     * incorrect face sampling without testing against a random ensemble mean.
     * The tolerance allows accumulated quadrature/MPI summation roundoff.
     */
    double h = L0/(1 << minlevel), reference = 0.;
    for (int m = 0; m < 9; ++m)
      for (int d = 0; d < 3; ++d)
        reference += 2.*(sq(forcing.transverse[m][d]) +
                         sq(forcing.transverse[m][3+d]))*
                     sq(cos(forcing.k[m][d]*h/2.));
    validation_require(isfinite(reference) &&
                       fabs(force2/volume - reference) <= 1e-8*reference,
                       "uniform-grid forcing violates filtered Parseval identity");
    validation_require(isfinite(mean_fx) && isfinite(mean_fy) && isfinite(mean_fz) &&
                       sqrt(sq(mean_fx/volume) + sq(mean_fy/volume) +
                            sq(mean_fz/volume)) <= 1e-8*sqrt(reference),
                       "uniform-grid forcing has a nonzero spatial mean");
    /** Only for this mode band: every nonzero k_d has the same magnitude,
     * so the face-difference symbol 2 sin(k_d h/2)/h is proportional to k_d.
     * Transversality therefore also implies discrete divergence-free forcing.
     * This roundoff check must not be applied across AMR interfaces or assumed
     * valid for an arbitrary wider forcing band.
     */
    validation_require(sqrt(div2/volume)*h <= 1e-8*sqrt(reference),
                       "uniform-grid forcing has nonzero discrete divergence");
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
