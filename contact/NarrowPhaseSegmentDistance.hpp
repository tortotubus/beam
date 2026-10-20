#pragma once

#include "elff/contact/NarrowPhase.hpp"
#include "elff/config/config.hpp"
#include <array>

namespace ELFF::Contact {

/** @brief Closed line segment in world coordinates; coincident endpoints are allowed. */
struct LineSegment {
  std::array<real_t, 3> start;
  std::array<real_t, 3> end;
};

/** @brief Closest points and interpolation parameters for two line segments. */
struct SegmentDistanceResult {
  real_t distance_squared; /**< @brief Squared centerline separation. */
  std::array<real_t, 3> point_a; /**< @brief Closest point on the first segment. */
  std::array<real_t, 3> point_b; /**< @brief Closest point on the second segment. */
  real_t parameter_a; /**< @brief First segment parameter in \f([0,1]\f). */
  real_t parameter_b; /**< @brief Second segment parameter in \f([0,1]\f). */
};

/**
 * @brief Closest-distance and clearance queries for piecewise-linear centerlines.
 *
 * For each segment pair, minimize
 * \f[
 *   \|\mathbf{a}_0+s(\mathbf{a}_1-\mathbf{a}_0)
 *     -\mathbf{b}_0-t(\mathbf{b}_1-\mathbf{b}_0)\|^2,
 *   \qquad s,t\in[0,1].
 * \f]
 * The solver checks the four boundary minima and the interior stationary point.
 * Parallel segments and zero-length segments are supported. No contact normal
 * is inferred at zero distance, where a unique normal need not exist.
 * Apply this query to segment pairs supplied by the broad phase; it does not
 * own curves, update the grid or compute contact forces.
 */
class NarrowPhaseSegmentDistance : public NarrowPhase {
public:
  /**
   * @brief Compute closest points, their parameters and squared separation.
   * @param a First finite segment.
   * @param b Second finite segment.
   * @return One minimizing pair; the minimizer need not be unique.
   * Invalid coordinates or an unrepresentable squared distance use ELFF_ABORT.
   */
  SegmentDistanceResult closest_points(const LineSegment& a,
                                       const LineSegment& b) const;

  /**
   * @brief Return whether centerline distance is at most the given clearance.
   * @param a First finite segment.
   * @param b Second finite segment.
   * @param clearance Finite nonnegative centerline exclusion distance.
   * For constant-radius fibres, pass the sum of both radii plus the required
   * surface gap. Equality is treated conservatively as a clearance violation.
   * Broad-phase bounds must include the corresponding search margin; this
   * method cannot recover pairs omitted by candidate generation.
   */
  bool within_clearance(const LineSegment& a, const LineSegment& b,
                        real_t clearance) const;
};

} // namespace ELFF::Contact
