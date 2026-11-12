#include "elff/contact/SphereGeometry.hpp"
#include "elff/general/error.hpp"
#include <cmath>

namespace ELFF::Contact {
namespace {
void validate(const SphereGeometry& sphere)
{
  ELFF_VERIFY(std::isfinite(sphere.radius) && sphere.radius > 0,
              "Sphere radius must be finite and positive.");
  for (auto x : sphere.center)
    ELFF_VERIFY(std::isfinite(x), "Sphere center must be finite.");
}
}

PlaneGeometry make_plane(const Vec3& point, const Vec3& normal)
{
  for (size_t d = 0; d < 3; ++d)
    ELFF_VERIFY(std::isfinite(point[d]) && std::isfinite(normal[d]),
                "Plane data must be finite.");
  const real_t length = std::hypot(normal[0], normal[1], normal[2]);
  ELFF_VERIFY(std::isfinite(length) && length > 0, "Plane normal must be nonzero.");
  PlaneGeometry plane{point, normal};
  for (auto& x : plane.normal) x /= length;
  return plane;
}

AABB sphere_bounds(const SphereGeometry& sphere)
{
  validate(sphere);
  AABB box;
  for (size_t d = 0; d < 3; ++d) {
    box.min[d] = sphere.center[d] - sphere.radius;
    box.max[d] = sphere.center[d] + sphere.radius;
    ELFF_VERIFY(std::isfinite(box.min[d]) && std::isfinite(box.max[d]),
                "Sphere bounds overflow.");
  }
  return box;
}

ContactPoint sphere_sphere(const SphereGeometry& a, const SphereGeometry& b)
{
  validate(a);
  validate(b);
  Vec3 normal;
  for (size_t d = 0; d < 3; ++d) normal[d] = a.center[d] - b.center[d];
  const real_t distance = std::hypot(normal[0], normal[1], normal[2]);
  ELFF_VERIFY(std::isfinite(distance) && distance > 0,
              "Coincident sphere centers or overflowing distance.");
  ContactPoint result{distance - a.radius - b.radius, {}, {}};
  ELFF_VERIFY(std::isfinite(result.gap), "Sphere gap overflow.");
  for (size_t d = 0; d < 3; ++d) {
    result.normal[d] = normal[d] / distance;
    // Midpoint between facing surface witnesses, without summing large centers.
    result.point[d] = a.center[d]
      - (a.radius + result.gap / 2) * result.normal[d];
    ELFF_VERIFY(std::isfinite(result.point[d]), "Contact point overflow.");
  }
  return result;
}

ContactPoint sphere_plane(const SphereGeometry& sphere, const PlaneGeometry& input)
{
  validate(sphere);
  const auto plane = make_plane(input.point, input.normal);
  real_t distance = 0;
  for (size_t d = 0; d < 3; ++d)
    distance += (sphere.center[d] - plane.point[d]) * plane.normal[d];
  ContactPoint result{distance - sphere.radius, plane.normal, {}};
  ELFF_VERIFY(std::isfinite(result.gap), "Sphere-plane gap overflow.");
  for (size_t d = 0; d < 3; ++d) {
    result.point[d] = sphere.center[d]
      - (sphere.radius + result.gap / 2) * plane.normal[d];
    ELFF_VERIFY(std::isfinite(result.point[d]), "Contact point overflow.");
  }
  return result;
}

} // namespace ELFF::Contact
