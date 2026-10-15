#pragma once

#include <cmath>
#include <cstddef>
#include <stdexcept>

#include <elff/models/beam/EulerBeamMesh.hpp>

namespace ELFF {

using namespace Models;

/**
 * @brief Manufactured dynamic solution for the shifted/pinned circular-arc
 * benchmark.
 *
 * This version pins the left endpoint at the origin:
 *
 *     \f(x(s,t) = 1/\phi(t) * \left[\begin{array}{c}\cos{(s \phi(t))} - 1
 * \\ \sin{(s \phi(t))} \\ 0\end{array}\right]\f),
 *
 * with \f(\phi(t) = \exp{(t)},\: s \in [0, l]\f).
 *
 * The tangent/slope is unchanged from the unshifted circular arc:
 *
 *     \f(x_s(s,t) = \left[\begin{array}{c}-\sin{(s \phi(t))} \\ \cos{(s
 * \phi(t))} \\ 0 \end{array}\right]\f).
 *
 * The velocity is
 *
 *     \f(x_t(s,t) = \left[ \begin{array}{c}(1 - \cos{(s \phi)})/\phi - s
 * \sin{(s \phi)}\\ -sin{(s \phi)}/\phi + s \cos{(s \phi)}\\ 0
 * \end{array}\right]\f),
 *
 * and the acceleration is
 *
 *     \f(x_{tt}(s,t) = \left[ \begin{array}{c}s \sin{(s \phi)} -
 * (1 - \cos{(s \phi)})/\phi - s^2 \phi \cos{(s \phi)}\\
 * \sin{(s \phi)}/\phi - s \cos{(s \phi)} - s^2 \phi \sin{(s \phi)}\\
 * 0 \end{array}\right]\f),
 *
 * for \f(\phi(t) = \exp{(t)}\right.
 *
 * @param N The number of nodes to compute; must be at least 2.
 * @param t The time for which to compute the solution.
 * @param length The physical length of the beam; default is \f(\pi/2\f).
 */
inline EulerBeamMesh
ManufacturedDynamicResult1(int N, double t, double length = M_PI_2)
{
  if (N < 2) {
    throw std::invalid_argument(
      "ManufacturedDynamicResult1 requires N >= 2 nodes.");
  }

  EulerBeamMesh mesh(N, length);

  auto& position = mesh.get_centerline();
  auto& velocity = mesh.get_centerline_velocity();
  auto& acceleration = mesh.get_centerline_acceleration();
  auto& slope = mesh.get_slope();
  auto& slope_velocity = mesh.get_slope_velocity();
  auto& slope_acceleration = mesh.get_slope_acceleration();

  const double phi = std::exp(t);
  const double ds = mesh.get_ds();

  for (int i = 0; i < N; ++i) {
    const std::size_t idx = static_cast<std::size_t>(i);

    const double s = static_cast<double>(i) * ds;
    const double theta = s * phi;

    const double c = std::cos(theta);
    const double sn = std::sin(theta);

    // Shifted/pinned position:
    //
    // x = (1/phi) * [cos(theta) - 1, sin(theta), 0]
    //
    // This gives x(0,t) = [0,0,0].

    position[idx][0] = (c - 1.0) / phi;
    position[idx][1] = sn / phi;
    position[idx][2] = 0.0;

    // Tangent / slope:
    //
    // x_s = [-sin(theta), cos(theta), 0]
    //
    // Translation does not change x_s.

    slope[idx][0] = -sn;
    slope[idx][1] = c;
    slope[idx][2] = 0.0;

    // Time derivative of the tangent / slope:
    //
    // x_st = [-s phi cos(theta), -s phi sin(theta), 0]

    slope_velocity[idx][0] = -s * phi * c;
    slope_velocity[idx][1] = -s * phi * sn;
    slope_velocity[idx][2] = 0.0;

    // Second time derivative of the tangent / slope:
    //
    // x_stt = [-s phi cos(theta) + s^2 phi^2 sin(theta),
    //          -s phi sin(theta) - s^2 phi^2 cos(theta),
    //          0]

    slope_acceleration[idx][0] = -s * phi * c + s * s * phi * phi * sn;
    slope_acceleration[idx][1] = -s * phi * sn - s * s * phi * phi * c;
    slope_acceleration[idx][2] = 0.0;

    // Velocity for the shifted solution:
    //
    // x_t = original_x_t - original_x_t(0,t)
    //
    // original_x_t = -(1/phi) e_r + s e_theta
    // original_x_t(0,t) = [-1/phi, 0, 0]

    velocity[idx][0] = (1.0 - c) / phi - s * sn;
    velocity[idx][1] = -sn / phi + s * c;
    velocity[idx][2] = 0.0;

    // Acceleration for the shifted solution:
    //
    // x_tt = [s sin(theta) - (1 - cos(theta))/phi - s^2 phi cos(theta),
    //         sin(theta)/phi - s cos(theta) - s^2 phi sin(theta),
    //         0]

    acceleration[idx][0] =
      s * sn - (1.0 - c) / phi - s * s * phi * c;
    acceleration[idx][1] = sn / phi - s * c - s * s * phi * sn;
    acceleration[idx][2] = 0.0;
  }

  return mesh;
}

/**
 * @brief Exact time derivative of the manufactured slope \f(x_s(s,t)\f).
 *
 * For
 *
 *   \f(x_s(s,t) = [-\sin(s \phi), \cos(s \phi), 0]\f)
 *
 * with \f(\phi(t) = \exp(t)\f), the time derivative is
 *
 *   \f(x_{st}(s,t) = [-s \phi \cos(s \phi), -s \phi \sin(s \phi), 0]\f).
 */
inline std::array<real_t, 3>
ManufacturedDynamicResult1SlopeVelocity(double s, double t)
{
  const double phi = std::exp(t);
  const double theta = s * phi;
  const double c = std::cos(theta);
  const double sn = std::sin(theta);

  return {-s * phi * c, -s * phi * sn, 0.0};
}

/**
 * @brief Exact second time derivative of the manufactured slope
 * \f(x_s(s,t)\f).
 *
 * For the same circular-arc benchmark, this is
 *
 *   \f(x_{stt}(s,t) = [-s \phi \cos(s \phi) + s^2 \phi^2 \sin(s \phi),
 *                      -s \phi \sin(s \phi) - s^2 \phi^2 \cos(s \phi), 0]\f).
 */
inline std::array<real_t, 3>
ManufacturedDynamicResult1SlopeAcceleration(double s, double t)
{
  const double phi = std::exp(t);
  const double theta = s * phi;
  const double c = std::cos(theta);
  const double sn = std::sin(theta);
  const double sphi = s * phi;
  const double ssphi2 = s * s * phi * phi;

  return {-sphi * c + ssphi2 * sn, -sphi * sn - ssphi2 * c, 0.0};
}

/**
 * @brief Manufactured linear forcing for ManufacturedDynamicResult1.
 *
 * This is the distributed forcing obtained by inserting the manufactured
 * solution into the linear beam equation
 *
 *   f_lin = mu * x_tt + EI * x_ssss.
 *
 * @param N Number of nodes where the forcing is sampled.
 * @param t Time at which the forcing is evaluated.
 * @param EI Bending stiffness.
 * @param mu Mass density per unit length.
 * @param length Physical beam length; default is pi/2.
 */
inline std::vector<std::array<real_t, 3>>
ManufacturedDynamicResult1LinearLoad(int    N,
                                     double t,
                                     double EI,
                                     double mu,
                                     double length = M_PI_2)
{
  if (N < 2) {
    throw std::invalid_argument(
      "ManufacturedDynamicResult1LinearLoad requires N >= 2 nodes.");
  }

  EulerBeamMesh mesh(N, length);
  const double  ds = mesh.get_ds();
  const double  phi = std::exp(t);

  std::vector<std::array<real_t, 3>> load(static_cast<std::size_t>(N),
                                          {0.0, 0.0, 0.0});

  for (int i = 0; i < N; ++i) {
    const std::size_t idx = static_cast<std::size_t>(i);
    const double      s = static_cast<double>(i) * ds;
    const double      theta = s * phi;
    const double      c = std::cos(theta);
    const double      sn = std::sin(theta);

    const double x_tt =
      s * sn - (1.0 - c) / phi - s * s * phi * c;
    const double y_tt = sn / phi - s * c - s * s * phi * sn;

    const double x_ssss = phi * phi * phi * c;
    const double y_ssss = phi * phi * phi * sn;

    load[idx][0] = mu * x_tt + EI * x_ssss;
    load[idx][1] = mu * y_tt + EI * y_ssss;
    load[idx][2] = 0.0;
  }

  return load;
}

/**
 * @brief Manufactured nonlinear forcing for ManufacturedDynamicResult1.
 *
 * This augments the linear manufactured force by the exact tangent-multiplier
 * contribution used in the nonlinear/inextensible benchmark:
 *
 *   f_nl = mu * x_tt + EI * x_ssss + (lambda * x_s)_s.
 *
 * For the circular-arc manufactured solution and constant lambda this becomes
 *
 *   f_nl = f_lin + lambda * x_ss.
 *
 * @param N Number of nodes where the forcing is sampled.
 * @param t Time at which the forcing is evaluated.
 * @param EI Bending stiffness.
 * @param mu Mass density per unit length.
 * @param lambda_value Constant manufactured tangent multiplier value.
 * @param length Physical beam length; default is pi/2.
 */
inline std::vector<std::array<real_t, 3>>
ManufacturedDynamicResult1NonlinearLoad(int    N,
                                        double t,
                                        double EI,
                                        double mu,
                                        double lambda_value = 1.0,
                                        double length = M_PI_2)
{
  auto load =
    ManufacturedDynamicResult1LinearLoad(N, t, EI, mu, length);

  EulerBeamMesh mesh(N, length);
  const double  ds = mesh.get_ds();
  const double  phi = std::exp(t);

  for (int i = 0; i < N; ++i) {
    const std::size_t idx = static_cast<std::size_t>(i);
    const double      s = static_cast<double>(i) * ds;
    const double      theta = s * phi;
    const double      c = std::cos(theta);
    const double      sn = std::sin(theta);

    const double x_ss = -phi * c;
    const double y_ss = -phi * sn;

    load[idx][0] += lambda_value * x_ss;
    load[idx][1] += lambda_value * y_ss;
  }

  return load;
}

} // namespace ELFF
