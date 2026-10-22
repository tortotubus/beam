#include <gtest/gtest.h>
#include "basilisk/templates/isotropic_forcing.h"

TEST(EswaranPope, ModesProjectionPeriodicityAndDeterminism) {
  EPForcing a, b;
  ASSERT_EQ(ep_init(&a, 2., 0.7, 0.3, 42, 1), 0);
  ASSERT_EQ(ep_init(&b, 2., 0.7, 0.3, 42, 1), 0);
  for (double dt : {0.01, 0.03, 1e-10, 0.2}) {
    ASSERT_EQ(ep_advance(&a, dt), 0);
    ASSERT_EQ(ep_advance(&b, dt), 0);
    for (int m = 0; m < 9; ++m) {
      int norm = 0;
      for (int d = 0; d < 3; ++d) norm += ep_wave_numbers[m][d]*ep_wave_numbers[m][d];
      EXPECT_GE(norm, 1); EXPECT_LE(norm, 2);
      for (int part = 0; part < 2; ++part) {
        double dot = 0.;
        for (int d = 0; d < 3; ++d) {
          dot += a.k[m][d]*a.transverse[m][3*part+d];
          EXPECT_EQ(a.b[m][3*part+d], b.b[m][3*part+d]);
        }
        EXPECT_NEAR(dot, 0., 1e-14);
      }
    }
    for (int d = 0; d < 3; ++d) {
      double f = ep_component(&a, .17, .28, .39, d);
      EXPECT_NEAR(f, ep_component(&a, 2.17, .28, .39, d), 1e-14);
      EXPECT_NEAR(f, ep_component(&a, .17, 2.28, .39, d), 1e-14);
      EXPECT_NEAR(f, ep_component(&a, .17, .28, 2.39, d), 1e-14);
    }
  }
  auto before = a.counter;
  EXPECT_EQ(ep_advance(&a, 0.), 0);
  EXPECT_EQ(a.counter, before);
  EXPECT_NE(ep_advance(&a, -1.), 0);
  EXPECT_EQ(a.counter, before);
}

TEST(EswaranPope, StationaryVarianceAndExactTransitionStatistics) {
  double sum = 0., square = 0., cross = 0., next_square = 0., physical_square = 0.;
  const int samples = 6000;
  const double sigma = .3, dt = .17, time = .7;
  for (int seed = 0; seed < samples; ++seed) {
    EPForcing state;
    ASSERT_EQ(ep_init(&state, 2., time, sigma, seed, 1), 0);
    for (int m = 0; m < 9; ++m)
      for (int d = 0; d < 6; ++d)
        physical_square += 2.*state.transverse[m][d]*state.transverse[m][d];
    auto previous = state;
    ASSERT_EQ(ep_advance(&state, dt), 0);
    for (int m = 0; m < 9; ++m)
      for (int d = 0; d < 6; ++d) {
        double x = previous.b[m][d], y = state.b[m][d];
        sum += x; square += x*x; cross += x*y; next_square += y*y;
      }
  }
  double count = samples*54.;
  EXPECT_NEAR(sum/count, 0., .002);
  EXPECT_NEAR(square/count, sigma*sigma, .001);
  EXPECT_NEAR(next_square/count, sigma*sigma, .001);
  EXPECT_NEAR(cross/square, exp(-dt/time), .004);
  EXPECT_NEAR(physical_square/samples, 72.*sigma*sigma, .05);
}

TEST(EswaranPope, ZeroAmplitudeAndInvalidControls) {
  EPForcing state;
  EXPECT_NE(ep_init(&state, 0., 1., 1., 1, 1), 0);
  EXPECT_NE(ep_init(&state, 1., 0., 1., 1, 1), 0);
  EXPECT_NE(ep_init(&state, 1., 1., -1., 1, 1), 0);
  ASSERT_EQ(ep_init(&state, 1., 1., 0., 1, 1), 0);
  ep_advance(&state, .1);
  for (int d = 0; d < 3; ++d) EXPECT_EQ(ep_component(&state, .1, .2, .3, d), 0.);
}
