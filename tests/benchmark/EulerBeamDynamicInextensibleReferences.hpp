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
 *     \f(x_{tt}(s,t) = \left[ \begin{array}{c}2 s \sin{(s \phi)} -
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
    // x_tt = [2 s sin(theta) - (1 - cos(theta))/phi - s^2 phi cos(theta),
    //         sin(theta)/phi - s cos(theta) - s^2 phi sin(theta),
    //         0]

    acceleration[idx][0] =
      2.0 * s * sn - (1.0 - c) / phi - s * s * phi * c;
    acceleration[idx][1] = sn / phi - s * c - s * s * phi * sn;
    acceleration[idx][2] = 0.0;
  }

  return mesh;
}

} // namespace ELFF
