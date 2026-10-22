#ifndef BASILISK_ESWARAN_POPE_H
#define BASILISK_ESWARAN_POPE_H
#include <math.h>
#include <stdint.h>
#include <string.h>

/** @brief Nine independent conjugate pairs in the band 0 < |n|^2 <= 2. */
static const int ep_wave_numbers[9][3] = {
  {1,0,0}, {0,1,0}, {0,0,1}, {1,1,0}, {1,-1,0},
  {1,0,1}, {1,0,-1}, {0,1,1}, {0,1,-1}
};

/** @brief OU state. b contains unprojected real then imaginary components.
 * sigma is the stationary standard deviation of each real component, not
 * physical-space RMS acceleration. The expected total RMS is 6 sqrt(2) sigma.
 */
typedef struct {
  double length, correlation_time, sigma;
  uint64_t seed, counter;
  double k[9][3], b[9][6], transverse[9][6];
} EPForcing;

/** @brief SplitMix64 mixer; unsigned overflow is intentional. */
static uint64_t ep_mix (uint64_t x) {
  x += UINT64_C(0x9e3779b97f4a7c15);
  x = (x ^ (x >> 30))*UINT64_C(0xbf58476d1ce4e5b9);
  x = (x ^ (x >> 27))*UINT64_C(0x94d049bb133111eb);
  return x ^ (x >> 31);
}

/** @brief Counter-indexed Box-Muller normal, independent of grid traversal. */
static double ep_normal (EPForcing *s) {
  uint64_t a = ep_mix(s->seed ^ ep_mix(s->counter++));
  uint64_t b = ep_mix(s->seed ^ ep_mix(s->counter++));
  double u = ((double)(a >> 11) + 1.)/9007199254740993.;
  double v = (double)(b >> 11)*0x1.0p-53;
  return sqrt(-2.*log(u))*cos(6.2831853071795864769*v);
}

/** @brief Project both complex components without amplitude normalization. */
static void ep_project (EPForcing *s) {
  for (int m = 0; m < 9; ++m) {
    double norm = 0.;
    for (int d = 0; d < 3; ++d) norm += ep_wave_numbers[m][d]*ep_wave_numbers[m][d];
    for (int part = 0; part < 2; ++part) {
      double dot = 0.;
      for (int d = 0; d < 3; ++d) dot += ep_wave_numbers[m][d]*s->b[m][3*part+d];
      for (int d = 0; d < 3; ++d)
        s->transverse[m][3*part+d] = s->b[m][3*part+d] - ep_wave_numbers[m][d]*dot/norm;
    }
  }
}

/** @brief Initialize with stationary Gaussian coefficients or a zero transient. */
static int ep_init (EPForcing *s, double length, double time, double sigma,
                    uint64_t seed, int stationary) {
  if (!isfinite(length) || length <= 0. || !isfinite(time) || time <= 0. ||
      !isfinite(sigma) || sigma < 0.) return 1;
  memset(s, 0, sizeof(*s));
  s->length = length; s->correlation_time = time; s->sigma = sigma; s->seed = seed;
  for (int m = 0; m < 9; ++m) {
    for (int d = 0; d < 3; ++d)
      s->k[m][d] = 6.2831853071795864769/length*ep_wave_numbers[m][d];
    if (stationary)
      for (int d = 0; d < 6; ++d) s->b[m][d] = sigma*ep_normal(s);
  }
  ep_project(s);
  return 0;
}

/** @brief Exact variable-step OU transition; zero dt consumes no random draws. */
static int ep_advance (EPForcing *s, double dt) {
  if (!isfinite(dt) || dt < 0.) return 1;
  if (dt == 0.) return 0;
  double decay = exp(-dt/s->correlation_time);
  double noise = s->sigma*sqrt(-expm1(-2.*dt/s->correlation_time));
  for (int m = 0; m < 9; ++m)
    for (int d = 0; d < 6; ++d)
      s->b[m][d] = decay*s->b[m][d] + noise*ep_normal(s);
  ep_project(s);
  return 0;
}

/** @brief Evaluate one acceleration component at coordinates relative to origin. */
static double ep_component (const EPForcing *s, double x, double y, double z, int d) {
  double value = 0.;
  for (int m = 0; m < 9; ++m) {
    double phase = s->k[m][0]*x + s->k[m][1]*y + s->k[m][2]*z;
    value += 2.*(s->transverse[m][d]*cos(phase) - s->transverse[m][3+d]*sin(phase));
  }
  return value;
}

#endif
