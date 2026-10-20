#include "elff/contact/NarrowPhaseSegmentDistance.hpp"
#include "elff/general/error.hpp"
#include <algorithm>
#include <cmath>
#include <limits>

namespace ELFF::Contact {
namespace {
using Vector = std::array<long double, 3>;
Vector subtract(const Vector& a, const Vector& b) {
  return {a[0]-b[0], a[1]-b[1], a[2]-b[2]};
}
long double dot(const Vector& a, const Vector& b) {
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}
Vector cross(const Vector& a, const Vector& b) {
  return {a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]};
}
Vector checked(const std::array<real_t, 3>& p) {
  for (auto x : p)
    if (!std::isfinite(x)) ELFF_ABORT("Segment coordinates must be finite.");
  return {p[0], p[1], p[2]};
}
struct Minimum {
  long double squared = std::numeric_limits<long double>::infinity();
  long double s = 0, t = 0;
  Vector a, b;
};
Minimum minimum(const LineSegment& first, const LineSegment& second) {
  const Vector a0 = checked(first.start), a1 = checked(first.end);
  const Vector b0 = checked(second.start), b1 = checked(second.end);
  Vector u = subtract(a1, a0), v = subtract(b1, b0), r = subtract(b0, a0);
  // Common scaling avoids large products in the minimization. Wider arithmetic
  // and the cross-product determinant retain nearly parallel interior minima.
  long double scale = 0;
  for (int d = 0; d < 3; ++d)
    scale = std::max({scale, std::abs(u[d]), std::abs(v[d]), std::abs(r[d])});
  if (!std::isfinite(scale)) ELFF_ABORT("Segment coordinate range is too large.");
  if (scale == 0) scale = 1;
  for (int d = 0; d < 3; ++d) { u[d] /= scale; v[d] /= scale; r[d] /= scale; }
  Minimum best;
  auto consider = [&](long double s, long double t) {
    s = std::clamp(s, 0.0L, 1.0L);
    t = std::clamp(t, 0.0L, 1.0L);
    Vector delta;
    for (int d = 0; d < 3; ++d) delta[d] = s*u[d] - t*v[d] - r[d];
    const long double squared = dot(delta, delta);
    if (squared < best.squared) { best.squared = squared; best.s = s; best.t = t; }
  };
  const long double uu = dot(u,u), vv = dot(v,v);
  const long double ur = dot(u,r), vr = dot(v,r), uv = dot(u,v);
  consider(0, vv > 0 ? -vr/vv : 0);
  consider(1, vv > 0 ? (uv-vr)/vv : 0);
  consider(uu > 0 ? ur/uu : 0, 0);
  consider(uu > 0 ? (ur+uv)/uu : 0, 1);
  const Vector normal = cross(u,v);
  const long double determinant = dot(normal,normal);
  if (determinant > 0) {
    const long double s = dot(cross(r,v),normal)/determinant;
    const long double t = dot(cross(r,u),normal)/determinant;
    if (s >= 0 && s <= 1 && t >= 0 && t <= 1) consider(s,t);
  }
  // Interpolate endpoints directly to avoid overflowing endpoint differences.
  for (int d = 0; d < 3; ++d) {
    best.a[d] = (1-best.s)*a0[d] + best.s*a1[d];
    best.b[d] = (1-best.t)*b0[d] + best.t*b1[d];
  }
  best.squared = (best.squared * scale) * scale;
  return best;
}
} // namespace

SegmentDistanceResult NarrowPhaseSegmentDistance::closest_points(
    const LineSegment& a, const LineSegment& b) const {
  const auto result = minimum(a,b);
  if (!std::isfinite(result.squared) || result.squared > std::numeric_limits<real_t>::max())
    ELFF_ABORT("Squared segment distance exceeds the configured real range.");
  return {static_cast<real_t>(result.squared),
          {static_cast<real_t>(result.a[0]), static_cast<real_t>(result.a[1]), static_cast<real_t>(result.a[2])},
          {static_cast<real_t>(result.b[0]), static_cast<real_t>(result.b[1]), static_cast<real_t>(result.b[2])},
          static_cast<real_t>(result.s), static_cast<real_t>(result.t)};
}

bool NarrowPhaseSegmentDistance::within_clearance(const LineSegment& a,
    const LineSegment& b, real_t clearance) const {
  if (!std::isfinite(clearance) || clearance < 0)
    ELFF_ABORT("Segment clearance must be finite and nonnegative.");
  const auto result = minimum(a,b);
  return std::sqrt(result.squared) <= static_cast<long double>(clearance);
}
} // namespace ELFF::Contact
