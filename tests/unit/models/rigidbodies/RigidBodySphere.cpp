#include "elff/models/rigidbody/RigidBodySphere.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

namespace ELFF {
namespace Models {
namespace {

using Vec3 = RigidBody3DModel::Vec3;

static bool
has_gnuplot()
{
  return std::system("gnuplot --version > /dev/null 2>&1") == 0;
}

static bool
has_display()
{
  const char* display = std::getenv("DISPLAY");
  const char* wayland = std::getenv("WAYLAND_DISPLAY");
  return (display && display[0] != '\0') || (wayland && wayland[0] != '\0');
}

TEST(RigidBodySphereTest, FibonacciGeometryAndMeasures)
{
  const real_t R = 2.0;
  const size_t N = 200;
  RigidBodySphere sphere(R, N, 1.0);

  const auto& pts = sphere.mesh().reference_points();
  const auto& nrm = sphere.mesh().reference_normals();
  const auto& ds = sphere.mesh().ds();

  ASSERT_EQ(pts.size(), N);
  ASSERT_EQ(nrm.size(), N);
  ASSERT_EQ(ds.size(), N);

  const real_t expected_ds = 4. * std::acos(-1.) * R * R / static_cast<real_t>(N);
  const real_t expected_area = 4. * std::acos(-1.) * R * R;

  for (size_t i = 0; i < N; ++i) {
    const real_t pnorm =
      std::sqrt(pts[i][0] * pts[i][0] + pts[i][1] * pts[i][1] + pts[i][2] * pts[i][2]);
    const real_t nnorm =
      std::sqrt(nrm[i][0] * nrm[i][0] + nrm[i][1] * nrm[i][1] + nrm[i][2] * nrm[i][2]);

    EXPECT_NEAR(pnorm, R, 1e-10);
    EXPECT_NEAR(nnorm, 1.0, 1e-10);
    EXPECT_NEAR(ds[i], expected_ds, 1e-12);
  }

  EXPECT_NEAR(sphere.mesh().total_measure(), expected_area, 1e-10);
}

TEST(RigidBodySphereTest, DensityControlsMassInStepResponse)
{
  const real_t R = 1.0;
  const size_t N = 100;
  const real_t rho = 2.0;
  RigidBodySphere sphere(R, N, rho);

  // Uniform traction in +x gives F = (4 pi R^2, 0, 0).
  std::vector<Vec3> traction(N, {1.0, 0.0, 0.0});
  const real_t dt = 0.1;
  sphere.step(dt, traction);

  const real_t m = rho * (4. / 3.) * std::acos(-1.) * R * R * R;
  const real_t F = 4. * std::acos(-1.) * R * R;
  const real_t v_expected = dt * (F / m);

  EXPECT_NEAR(sphere.com_velocity()[0], v_expected, 1e-10);
  EXPECT_NEAR(sphere.com_velocity()[1], 0.0, 1e-12);
  EXPECT_NEAR(sphere.com_velocity()[2], 0.0, 1e-12);
}

TEST(RigidBodySphereTest, PlotMeshWithGnuplot)
{
  if (!has_gnuplot()) {
    GTEST_SKIP() << "gnuplot not available";
  }
  if (!has_display()) {
    GTEST_SKIP() << "no DISPLAY/WAYLAND_DISPLAY available for qt terminal";
  }

  RigidBodySphere sphere(1.0, 1000, 1.0);
  const auto& pts = sphere.mesh().reference_points();

  FILE* pipe = popen("gnuplot -persist", "w");
  ASSERT_NE(pipe, nullptr);

  std::fprintf(pipe, "set term qt 3\n");
  std::fprintf(pipe, "set title 'Fibonacci Sphere Mesh'\n");
  std::fprintf(pipe, "set xlabel 'x'\nset ylabel 'y'\nset zlabel 'z'\n");
  std::fprintf(pipe, "splot '-' with points pt 7 ps 0.4 notitle\n");
  for (const auto& p : pts) {
    std::fprintf(pipe, "%.16e %.16e %.16e\n", p[0], p[1], p[2]);
  }
  std::fprintf(pipe, "e\n");

  const int rc = pclose(pipe);
  EXPECT_EQ(rc, 0);
}

} // namespace
} // namespace Models
} // namespace ELFF
