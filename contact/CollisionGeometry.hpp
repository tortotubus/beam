#pragma once
#include "elff/contact/NarrowPhaseSegmentDistance.hpp"
#include "elff/contact/SphereGeometry.hpp"
#include <variant>

namespace ELFF::Contact {

/** @brief Segment swept by radius; degenerate segments are spheres. */
struct CapsuleGeometry {
  LineSegment axis;
  real_t radius;
};

using CollisionShape = std::variant<SphereGeometry, CapsuleGeometry>;
/**
 * @brief Analytic collision description in body coordinates, independent of
 * IBM markers.
 */
class CollisionGeometryProvider {
public:
  virtual ~CollisionGeometryProvider() = default;
  virtual CollisionShape collision_shape() const = 0;
};

AABB collision_bounds(const CollisionShape& shape);
/** @brief Closest-axis penalty contact, regularized at zero capsule-axis separation.
 * The hint orients ambiguous normals, for example opposite relative velocity.
 * This is not a general minimum-translation penetration solver. Coincident point
 * shapes remain unsupported. Identical static capsules require stable caller
 * ordering to assign force signs reproducibly.
 */
ContactPoint collision_pair(const CollisionShape& a, const CollisionShape& b,
                            const Vec3& normal_hint = {});
ContactPoint collision_plane(const CollisionShape& shape, const PlaneGeometry& plane);

} // namespace ELFF::Contact
