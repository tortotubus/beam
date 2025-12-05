#pragma once

#include <cmath>
#include <limits>
#include <stdexcept>

namespace beam {

struct BisshoppAndDrucker1945Result
{
  double k;     // elliptic modules in [1/sqrt(2), 1)
  double a;     // dimensionless load a = sqrt(PL^2/B)
  double delta; // vertical drop
  double A;     // horizontal retreat (L - x_tip)
};

/*
  For reference see https://doi.org/10.1090/qam/13360
*/
inline BisshoppAndDrucker1945Result
BisshoppAndDrucker1945(double length, double EI, double tip_load)
{
  /* Note in the original paper,
    L = beam length
    B = flexural rigidity
    P = vertical tip load
    \delta = tip vertical drop
    A = tip horizontal displacement
    \phi = slope angle
    s = arc-length
  */

  double P = tip_load;
  double L = length;
  double B = EI;

  BisshoppAndDrucker1945Result result{};

  // Compute the dimensionless load
  result.a = std::sqrt((P * L * L) / B);

  // Check if the load is zero
  if (result.a == 0.) {
    result.k = M_SQRT1_2;
    result.delta = 0.;
    result.A = 0.;
    return result;
  }

  // Helper function: a(k) = K(k) - F(theta1 | k), with theta1 =
  // asin(1/sqrt(2)*k);
  auto a_of_k = [](double k) -> double {
    const double kmin = M_SQRT1_2;   // 1/sqrt(2)
    const double kmax = 1.0 - 1e-15; // ~1
    k = k < kmin ? kmin : (k > kmax ? kmax : k);
    const double theta1 = std::asin(1. / (std::sqrt(2.) * k));
    const double K = std::comp_ellint_1(k);
    const double F1 = std::ellint_1(k, theta1);
    return K - F1;
  };

  const double k_lo0 = M_SQRT1_2;
  const double k_hi0 = 1.0 - 1e-12;
  double alo = a_of_k(k_lo0);
  double ahi = a_of_k(k_hi0);

  if (!(alo <= result.a && result.a < ahi)) {
    throw std::runtime_error("Failed to bracket target 'a'.");
  }

  // Bisection search parameters
  size_t max_it = 1000;
  double tol = 1e-10;

  // Initial interval
  double k_lo = k_lo0, k_hi = k_hi0;

  // Perform a bisection search on a(k) - a_target
  for (size_t it = 0; it < 1000; it++) {
    const double k_i = 0.5 * (k_lo + k_hi);
    const double a_i = a_of_k(k_i);

    if (std::abs(a_i - result.a) < tol) {
      result.k = k_i;
      break;
    }

    if (a_i < result.a) {
      k_lo = k_i;
    } else {
      k_hi = k_i;
    }

    if (std::abs(k_hi - k_lo) < tol) {
      result.k = 0.5 * (k_lo + k_hi);
      break;
    }

    if (it == max_it - 1) {
      result.k = 0.5 * (k_lo + k_hi);
    }
  }

  // Now that k has been computed, compute delta/L and (L-A)/L
  const double k = result.k;
  const double theta1 = std::asin(1. / (std::sqrt(2.) * k));
  const double K = std::comp_ellint_1(k);
  const double E = std::comp_ellint_2(k);
  const double F1 = std::ellint_1(k, theta1);
  const double E1 = std::ellint_2(k, theta1);

  const double a_val = K - F1;
  // const double delta_over_L = (E-E1)/a_val;
  const double delta_over_L = 1.0 - 2.0 * (E - E1) / a_val;
  const double x_over_L = (std::sqrt(2.) * std::sqrt(2. * k * k - 1.)) / a_val;

  result.delta = delta_over_L * L;
  result.A = (1. - x_over_L) * L;

  return result;
}
}