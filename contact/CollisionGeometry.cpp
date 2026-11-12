#include "elff/contact/CollisionGeometry.hpp"
#include "elff/general/error.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <tuple>

namespace ELFF::Contact {
namespace {
CapsuleGeometry swept(const CollisionShape& shape) {
  if (const auto* sphere = std::get_if<SphereGeometry>(&shape))
    return {{sphere->center, sphere->center}, sphere->radius};
  return std::get<CapsuleGeometry>(shape);
}
Vec3 subtract(const Vec3& a, const Vec3& b) { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }
real_t dot(const Vec3& a, const Vec3& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
} // namespace

AABB collision_bounds(const CollisionShape& shape) {
  const auto cap = swept(shape);
  const auto a = sphere_bounds({cap.axis.start, cap.radius});
  const auto b = sphere_bounds({cap.axis.end, cap.radius});
  AABB box;
  for (size_t d = 0; d < 3; ++d) {
    box.min[d] = std::min(a.min[d], b.min[d]);
    box.max[d] = std::max(a.max[d], b.max[d]);
  }
  return box;
}

ContactPoint collision_pair(const CollisionShape& a, const CollisionShape& b,
                            const Vec3& normal_hint) {
  collision_bounds(a);
  collision_bounds(b);
  const auto ca = swept(a), cb = swept(b);
  NarrowPhaseSegmentDistance query;
  auto result = query.closest_points(ca.axis, cb.axis);
  const auto u = subtract(ca.axis.end, ca.axis.start);
  const auto v = subtract(cb.axis.end, cb.axis.start);
  const real_t lu = std::hypot(u[0], u[1], u[2]);
  const real_t lv = std::hypot(v[0], v[1], v[2]);
  // Choose the midpoint of the minimizing interval for parallel overlapping axes.
  if (lu > 0 && lv > 0) {
    Vec3 axis, other;
    for (size_t d = 0; d < 3; ++d) {
      axis[d] = u[d] / lu;
      other[d] = v[d] / lv;
    }
    const Vec3 cross{axis[1] * other[2] - axis[2] * other[1],
                     axis[2] * other[0] - axis[0] * other[2],
                     axis[0] * other[1] - axis[1] * other[0]};
    if (std::hypot(cross[0], cross[1], cross[2]) <= 32 * std::numeric_limits<real_t>::epsilon()) {
      const real_t t0 = dot(subtract(cb.axis.start, ca.axis.start), axis);
      const real_t t1 = dot(subtract(cb.axis.end, ca.axis.start), axis);
      const real_t lower = std::max(real_t{0}, std::min(t0, t1));
      const real_t upper = std::min(lu, std::max(t0, t1));
      if (lower <= upper) {
        const real_t middle = lower + (upper - lower) / 2;
        for (size_t d = 0; d < 3; ++d)
          result.point_a[d] = ca.axis.start[d] + middle * axis[d];
        const real_t along =
            std::clamp(dot(subtract(result.point_a, cb.axis.start), other), real_t{0}, lv);
        for (size_t d = 0; d < 3; ++d)
          result.point_b[d] = cb.axis.start[d] + along * other[d];
        result.parameter_a = middle / lu;
        result.parameter_b = along / lv;
      }
    }
  }
  const auto separation = subtract(result.point_a, result.point_b);
  const real_t distance = std::hypot(separation[0], separation[1], separation[2]);
  const real_t tolerance =
      32 * std::numeric_limits<real_t>::epsilon() * std::max({lu, lv, ca.radius, cb.radius});
  ContactPoint contact;
  if (distance > tolerance || (lu == 0 && lv == 0)) {
    // Keep the ordinary witness normal whenever it is resolved.
    contact = sphere_sphere({result.point_a, ca.radius}, {result.point_b, cb.radius});
  } else {
    for (auto value : normal_hint)
      ELFF_VERIFY(std::isfinite(value), "Contact normal hint must be finite.");
    // Endpoint ordering is a parametrization, not an orientation of the body.
    auto direction = [](const CapsuleGeometry& capsule, real_t length) {
      Vec3 axis{};
      if (length > 0) {
        const bool reverse = capsule.axis.end < capsule.axis.start;
        for (size_t d = 0; d < 3; ++d)
          axis[d] = (reverse ? capsule.axis.start[d] - capsule.axis.end[d]
                             : capsule.axis.end[d] - capsule.axis.start[d]) /
                    length;
      }
      return axis;
    };
    const auto x = direction(ca, lu), y = direction(cb, lv);
    Vec3 normal{x[1] * y[2] - x[2] * y[1], x[2] * y[0] - x[0] * y[2], x[0] * y[1] - x[1] * y[0]};
    real_t norm = std::hypot(normal[0], normal[1], normal[2]);
    if (norm > 32 * std::numeric_limits<real_t>::epsilon()) {
      for (auto& value : normal)
        value /= norm;
      if (dot(normal, normal_hint) < 0)
        for (auto& value : normal)
          value = -value;
    } else {
      const auto& axis = lu > 0 ? x : y;
      const real_t hint_length = std::hypot(normal_hint[0], normal_hint[1], normal_hint[2]);
      Vec3 hint{};
      if (std::isfinite(hint_length) && hint_length > 0)
        for (size_t d = 0; d < 3; ++d)
          hint[d] = normal_hint[d] / hint_length;
      const real_t along = dot(hint, axis);
      for (size_t d = 0; d < 3; ++d)
        normal[d] = hint[d] - along * axis[d];
      norm = std::hypot(normal[0], normal[1], normal[2]);
      if (norm > 32 * std::numeric_limits<real_t>::epsilon()) {
        for (auto& value : normal)
          value /= norm;
      } else {
        // A deterministic transverse direction, signed by canonical geometry.
        size_t basis = 0;
        for (size_t d = 1; d < 3; ++d)
          if (std::abs(axis[d]) < std::abs(axis[basis]))
            basis = d;
        Vec3 e{};
        e[basis] = 1;
        normal = {axis[1] * e[2] - axis[2] * e[1], axis[2] * e[0] - axis[0] * e[2],
                  axis[0] * e[1] - axis[1] * e[0]};
        norm = std::hypot(normal[0], normal[1], normal[2]);
        auto key = [](const CapsuleGeometry& capsule) {
          return std::make_tuple(std::min(capsule.axis.start, capsule.axis.end),
                                 std::max(capsule.axis.start, capsule.axis.end), capsule.radius);
        };
        const real_t sign = key(cb) < key(ca) ? real_t{-1} : real_t{1};
        for (auto& value : normal)
          value *= sign / norm;
      }
    }
    contact.gap = distance - ca.radius - cb.radius;
    contact.normal = normal;
    contact.normal_fallback = true;
    for (size_t d = 0; d < 3; ++d) {
      const real_t middle = result.point_a[d] + (result.point_b[d] - result.point_a[d]) / 2;
      contact.point[d] = middle + (cb.radius - ca.radius) * normal[d] / 2;
      ELFF_VERIFY(std::isfinite(contact.point[d]), "Regularized contact point overflow.");
    }
    ELFF_VERIFY(std::isfinite(contact.gap), "Regularized gap overflow.");
  }
  contact.parameter_a = result.parameter_a;
  contact.parameter_b = result.parameter_b;
  return contact;
}

ContactPoint collision_plane(const CollisionShape& shape, const PlaneGeometry& input) {
  collision_bounds(shape);
  const auto cap = swept(shape);
  const auto plane = make_plane(input.point, input.normal);
  const real_t a = dot(subtract(cap.axis.start, plane.point), plane.normal);
  const real_t b = dot(subtract(cap.axis.end, plane.point), plane.normal);
  const auto delta = subtract(cap.axis.end, cap.axis.start);
  const real_t length = std::hypot(delta[0], delta[1], delta[2]);
  real_t parameter = a < b ? real_t{0} : real_t{1};
  if (std::abs(a - b) <= 32 * std::numeric_limits<real_t>::epsilon() * length)
    parameter = length > 0 ? real_t{0.5} : real_t{0};
  Vec3 point;
  for (size_t d = 0; d < 3; ++d)
    point[d] = cap.axis.start[d] + parameter * delta[d];
  auto contact = sphere_plane({point, cap.radius}, plane);
  contact.parameter_a = parameter;
  return contact;
}
} // namespace ELFF::Contact
