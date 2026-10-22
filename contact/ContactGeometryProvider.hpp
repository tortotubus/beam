#pragma once
#include "elff/config/config.hpp"
#include <array>
#include <cstddef>

namespace ELFF::Contact {
/** @brief Borrowed world-space polyline; consecutive points define segments.
 * Reacquire after initialization, advancement or restoration. No thickness or
 * ownership is implied. Repeated points represent valid degenerate segments.
 */
struct PolylineView {
  const std::array<real_t, 3> *points = nullptr;
  std::size_t count = 0;
};

/** @brief Optional geometric capability for clearance and future contact
 * queries. Beam providers expose nodal polylines, not exact curved-element
 * interpolation.
 */
class ContactGeometryProvider {
public:
  virtual ~ContactGeometryProvider() = default;
  /** @brief Borrow the current structural centerline in global coordinates. */
  virtual PolylineView contact_geometry() const = 0;
};
} // namespace ELFF::Contact
