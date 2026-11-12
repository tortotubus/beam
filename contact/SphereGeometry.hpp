#pragma once

#include "elff/contact/BroadPhase.hpp"
#include "elff/contact/Wrench.hpp"

namespace ELFF::Contact {

struct SphereGeometry {
  Vec3 center;
  real_t radius;
};

/** @brief Fixed infinite plane; normal points into the permitted half-space.
 * Reaction torque is reported about point. Construct through make_plane().
 */
struct PlaneGeometry {
  Vec3 point;
  Vec3 normal;
};

/** @brief Signed gap and common force application point.
 * Negative gap means penetration; normal is the compressive force direction on A.
 */
struct ContactPoint {
  real_t gap;
  Vec3 normal;
  Vec3 point;
  real_t parameter_a = 0; ///< Material coordinate on shape A, in [0,1].
  real_t parameter_b = 0; ///< Material coordinate on shape B, in [0,1].
  bool normal_fallback = false; ///< Normal chosen by the zero-distance regularization.
};

PlaneGeometry make_plane(const Vec3& point, const Vec3& normal);
AABB sphere_bounds(const SphereGeometry& sphere);
ContactPoint sphere_sphere(const SphereGeometry& a, const SphereGeometry& b);
ContactPoint sphere_plane(const SphereGeometry& sphere, const PlaneGeometry& plane);

} // namespace ELFF::Contact
