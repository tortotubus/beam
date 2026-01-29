/**
 * @file Shapes.hpp
 * @brief Finite element shape function library.
 *
 * This module provides compile-time templates for various shape functions
 * commonly used in finite element methods (FEM), including:
 * - Cubic-Hermite (C1-continuous) polynomials for Euler-Bernoulli beams
 * - Linear hat functions for standard FEM
 * - Quadratic Lagrange polynomials for higher-order elements
 *
 * All shape functions are defined on the reference element \f( [0,1] \f)
 * with template parameter \f( T \f) allowing automatic differentiation
 * or other numeric types.
 */

#pragma once
#include <array>
#include <cstddef>

namespace ELFF {

/**
 * @class CubicHermite
 * 
 * @brief Cubic-Hermite shape functions on the reference element \f( [0,1] \f).
 *
 * Cubic-Hermite polynomials are C1-continuous piecewise cubics defined by
 * nodal displacements (H1, H3) and nodal slopes/rotations (H2, H4).
 * These are the standard basis functions for Euler-Bernoulli beam theory.
 *
 * For an element with length \f( h \f), the interpolation is:
 * \f[ u(x) = H_1(\xi) u_1 + H_2(\xi, h) \theta_1 + H_3(\xi) u_2 + H_4(\xi, h) \theta_2 \f]
 *
 * where \f( \xi = x/h \in [0,1] \f) is the normalized reference coordinate.
 *
 * @tparam T Numeric type (e.g., double, float, or autodiff type)
 */
template<class T>
struct CubicHermite
{
  /**
   * @name Nodal Displacement Shape Functions
   * @brief First and second nodal displacement basis functions.
   * @{
   */

  /**
   * @brief Shape function \f( H_1(\xi) \f) for displacement at node 1.
   *
   * Returns 1 at \f( \xi = 0 \f) and 0 at \f( \xi = 1 \f).
   * First derivative vanishes at both ends.
   *
   * @param xi Reference coordinate in \f( [0,1] \f)
   * @return Value of \f( H_1 = 1 - 3\xi^2 + 2\xi^3 \f)
   */
  static T H1(T xi) noexcept
  {
    return 1.0 - 3.0 * xi * xi + 2.0 * xi * xi * xi;
  }

  /**
   * @brief Shape function \f( H_3(\xi) \f) for displacement at node 2.
   *
   * Returns 0 at \f( \xi = 0 \f) and 1 at \f( \xi = 1 \f).
   * First derivative vanishes at both ends.
   *
   * @param xi Reference coordinate in \f( [0,1] \f)
   * @return Value of \f( H_3 = 3\xi^2 - 2\xi^3 \f)
   */
  static T H3(T xi) noexcept { return 3.0 * xi * xi - 2.0 * xi * xi * xi; }

  /** @} */

  /**
   * @name Nodal Slope (Rotation) Shape Functions
   * @brief First and second nodal slope/rotation basis functions.
   * @{
   */

  /**
   * @brief Shape function \f( H_2(\xi, h) \f) for rotation at node 1.
   *
   * Returns 0 at \f( \xi = 0 \f) and 0 at \f( \xi = 1 \f).
   * First derivative equals 1 at \f( \xi = 0 \f) and 0 at \f( \xi = 1 \f).
   * Scaled by element length \f( h \f) for physical dimensions.
   *
   * @param xi Reference coordinate in \f( [0,1] \f)
   * @param h  Element length (scaling factor)
   * @return Value of \f( H_2 = h(\xi - 2\xi^2 + \xi^3) \f)
   */
  static T H2(T xi, T h) noexcept
  {
    return h * (xi - 2.0 * xi * xi + xi * xi * xi);
  }

  /**
   * @brief Shape function \f( H_4(\xi, h) \f) for rotation at node 2.
   *
   * Returns 0 at \f( \xi = 0 \f) and 0 at \f( \xi = 1 \f).
   * First derivative equals 0 at \f( \xi = 0 \f) and 1 at \f( \xi = 1 \f).
   * Scaled by element length \f( h \f) for physical dimensions.
   *
   * @param xi Reference coordinate in \f( [0,1] \f)
   * @param h  Element length (scaling factor)
   * @return Value of \f( H_4 = h(\xi^3 - \xi^2) \f)
   */
  static T H4(T xi, T h) noexcept { return h * (xi * xi * xi - xi * xi); }

  /** @} */

  /**
   * @name First Derivatives (w.r.t. Reference Coordinate)
   * @brief Derivatives \f( \frac{d}{d\xi} \f) of shape functions.
   * @{
   */

  /**
   * @brief First derivative of \f( H_1 \f) w.r.t. \f( \xi \f).
   * @param xi Reference coordinate in \f( [0,1] \f)
   * @return \f( \frac{dH_1}{d\xi} = -6\xi + 6\xi^2 \f)
   */
  static T dH1_dxi(T xi) noexcept { return -6.0 * xi + 6.0 * xi * xi; }

  /**
   * @brief First derivative of \f( H_2 \f) w.r.t. \f( \xi \f).
   * @param xi Reference coordinate in \f( [0,1] \f)
   * @param h  Element length (scaling factor)
   * @return \f( \frac{dH_2}{d\xi} = h(1 - 4\xi + 3\xi^2) \f)
   */
  static T dH2_dxi(T xi, T h) noexcept
  {
    return h * (1.0 - 4.0 * xi + 3.0 * xi * xi);
  }

  /**
   * @brief First derivative of \f( H_3 \f) w.r.t. \f( \xi \f).
   * @param xi Reference coordinate in \f( [0,1] \f)
   * @return \f( \frac{dH_3}{d\xi} = 6\xi - 6\xi^2 \f)
   */
  static T dH3_dxi(T xi) noexcept { return 6.0 * xi - 6.0 * xi * xi; }

  /**
   * @brief First derivative of \f( H_4 \f) w.r.t. \f( \xi \f).
   * @param xi Reference coordinate in \f( [0,1] \f)
   * @param h  Element length (scaling factor)
   * @return \f( \frac{dH_4}{d\xi} = h(3\xi^2 - 2\xi) \f)
   */
  static T dH4_dxi(T xi, T h) noexcept
  {
    return h * (3.0 * xi * xi - 2.0 * xi);
  }

  /** @} */

  /**
   * @name Second Derivatives (w.r.t. Reference Coordinate)
   * @brief Second derivatives \f( \frac{d^2}{d\xi^2} \f) of shape functions.
   * @{
   */

  /**
   * @brief Second derivative of \f( H_1 \f) w.r.t. \f( \xi \f).
   * @param xi Reference coordinate in \f( [0,1] \f)
   * @return \f( \frac{d^2H_1}{d\xi^2} = -6 + 12\xi \f)
   */
  static T d2H1_dxi2(T xi) noexcept { return -6.0 + 12.0 * xi; }

  /**
   * @brief Second derivative of \f( H_2 \f) w.r.t. \f( \xi \f).
   * @param xi Reference coordinate in \f( [0,1] \f)
   * @param h  Element length (scaling factor)
   * @return \f( \frac{d^2H_2}{d\xi^2} = h(-4 + 6\xi) \f)
   */
  static T d2H2_dxi2(T xi, T h) noexcept { return h * (-4.0 + 6.0 * xi); }

  /**
   * @brief Second derivative of \f( H_3 \f) w.r.t. \f( \xi \f).
   * @param xi Reference coordinate in \f( [0,1] \f)
   * @return \f( \frac{d^2H_3}{d\xi^2} = 6 - 12\xi \f)
   */
  static T d2H3_dxi2(T xi) noexcept { return 6.0 - 12.0 * xi; }

  /**
   * @brief Second derivative of \f( H_4 \f) w.r.t. \f( \xi \f).
   * @param xi Reference coordinate in \f( [0,1] \f)
   * @param h  Element length (scaling factor)
   * @return \f( \frac{d^2H_4}{d\xi^2} = h(6\xi - 2) \f)
   */
  static T d2H4_dxi2(T xi, T h) noexcept { return h * (6.0 * xi - 2.0); }

  /** @} */

  /**
   * @name Utility Functions
   * @brief Convenience functions that return multiple values at once.
   * @{
   */

  /**
   * @brief Evaluate all four shape functions at once.
   *
   * @param xi Reference coordinate in \f( [0,1] \f)
   * @param h  Element length (scaling factor)
   * @return Array \f( [H_1, H_2, H_3, H_4] \f)
   */
  static std::array<T, 4> values(T xi, T h) noexcept
  {
    return { H1(xi), H2(xi, h), H3(xi), H4(xi, h) };
  }

  /**
   * @brief Evaluate all four first derivatives (w.r.t. physical \f( x \f)).
   *
   * Computes \f( \frac{d}{dx} = \frac{1}{h} \frac{d}{d\xi} \f).
   *
   * @param xi Reference coordinate in \f( [0,1] \f)
   * @param h  Element length
   * @return Array of \f( [\frac{dH_1}{dx}, \frac{dH_2}{dx}, \frac{dH_3}{dx}, \frac{dH_4}{dx}] \f)
   */
  static std::array<T, 4> derivs(T xi, T h) noexcept
  {
    return {
      dH1_dxi(xi) / h, dH2_dxi(xi, h) / h, dH3_dxi(xi) / h, dH4_dxi(xi, h) / h
    };
  }

  /**
   * @brief Evaluate all four second derivatives (w.r.t. physical \f( x \f)).
   *
   * Computes \f( \frac{d^2}{dx^2} = \frac{1}{h^2} \frac{d^2}{d\xi^2} \f).
   *
   * @param xi Reference coordinate in \f( [0,1] \f)
   * @param h  Element length
   * @return Array of \f( [\frac{d^2H_1}{dx^2}, \frac{d^2H_2}{dx^2}, \frac{d^2H_3}{dx^2}, \frac{d^2H_4}{dx^2}] \f)
   */
  static std::array<T, 4> second_derivs(T xi, T h) noexcept
  {
    return { d2H1_dxi2(xi) / (h * h),
             d2H2_dxi2(xi, h) / (h * h),
             d2H3_dxi2(xi) / (h * h),
             d2H4_dxi2(xi, h) / (h * h) };
  }

  /** @} */
};

/**
 * @class LinearShape
 * @brief Linear (hat) shape functions on the reference element \f( [0,1] \f).
 *
 * These are the simplest piecewise-linear basis functions for standard
 * finite element methods. Each element is defined by two nodal values.
 *
 * Interpolation: \f[ u(x) = M_1(\xi) u_1 + M_2(\xi) u_2 \f]
 * where \f( \xi = x/h \in [0,1] \f).
 *
 * @tparam T Numeric type (e.g., double, float)
 */
template<class T>
struct LinearShape
{
  /**
   * @brief Shape function \f( M_1(\xi) \f) for node 1.
   *
   * Returns 1 at \f( \xi = 0 \f) and 0 at \f( \xi = 1 \f).
   *
   * @param xi Reference coordinate in \f( [0,1] \f)
   * @return Value of \f( M_1 = 1 - \xi \f)
   */
  static T M1(T xi) noexcept { return 1.0 - xi; }

  /**
   * @brief Shape function \f( M_2(\xi) \f) for node 2.
   *
   * Returns 0 at \f( \xi = 0 \f) and 1 at \f( \xi = 1 \f).
   *
   * @param xi Reference coordinate in \f( [0,1] \f)
   * @return Value of \f( M_2 = \xi \f)
   */
  static T M2(T xi) noexcept { return xi; }

  /**
   * @brief Evaluate both shape functions at once.
   *
   * @param xi Reference coordinate in \f( [0,1] \f)
   * @return Array \f( [M_1, M_2] \f)
   */
  static std::array<T, 2> values(T xi) noexcept { return { M1(xi), M2(xi) }; }
};

/**
 * @class QuadraticLagrange
 * @brief Quadratic Lagrange shape functions on the reference element \f( [0,1] \f).
 *
 * Quadratic Lagrange polynomials are second-order elements with three nodes:
 * two at the ends and one at the midpoint. They provide higher accuracy
 * than linear elements for smooth solutions.
 *
 * Interpolation: \f[ u(x) = L_0(\xi) u_1 + L_m(\xi) u_m + L_1(\xi) u_2 \f]
 * where \f( \xi = x/h \in [0,1] \f), and the midpoint node is at \f( \xi = 0.5 \f).
 *
 * @tparam T Numeric type (e.g., double, float)
 */
template<class T>
struct QuadraticLagrange
{
  /**
   * @brief Shape function \f( L_0(\xi) \f) for node at \f( \xi = 0 \f).
   *
   * Returns 1 at \f( \xi = 0 \f), 0 at \f( \xi = 0.5 \f), and 0 at \f( \xi = 1 \f).
   *
   * @param xi Reference coordinate in \f( [0,1] \f)
   * @return Value of \f( L_0 = 2(\xi - 0.5)(\xi - 1) \f)
   */
  static T L0(T xi) noexcept { return 2.0 * (xi - 0.5) * (xi - 1.0); }

  /**
   * @brief Shape function \f( L_m(\xi) \f) for midpoint node at \f( \xi = 0.5 \f).
   *
   * Returns 0 at \f( \xi = 0 \f), 1 at \f( \xi = 0.5 \f), and 0 at \f( \xi = 1 \f).
   * This is the only shape function that is non-zero at the midpoint.
   *
   * @param xi Reference coordinate in \f( [0,1] \f)
   * @return Value of \f( L_m = 4\xi(1 - \xi) \f)
   */
  static T Lm(T xi) noexcept { return 4.0 * xi * (1.0 - xi); }

  /**
   * @brief Shape function \f( L_1(\xi) \f) for node at \f( \xi = 1 \f).
   *
   * Returns 0 at \f( \xi = 0 \f), 0 at \f( \xi = 0.5 \f), and 1 at \f( \xi = 1 \f).
   *
   * @param xi Reference coordinate in \f( [0,1] \f)
   * @return Value of \f( L_1 = 2\xi(\xi - 0.5) \f)
   */
  static T L1(T xi) noexcept { return 2.0 * xi * (xi - 0.5); }

  /**
   * @brief Evaluate all three shape functions at once.
   *
   * @param xi Reference coordinate in \f( [0,1] \f)
   * @return Array \f( [L_0, L_m, L_1] \f)
   */
  static std::array<T, 3> values(T xi) noexcept
  {
    return { L0(xi), Lm(xi), L1(xi) };
  }
};

}
