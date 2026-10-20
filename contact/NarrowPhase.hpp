#pragma once

namespace ELFF::Contact {

/**
 * @brief Computes contact geometry for a pair of collision shapes.
 *
 * Implementations should evaluate candidate pairs from BroadPhase and return
 * signed separation or penetration, contact normals and witness points.
 * Geometry backends, such as GJK with a penetration solver, should remain
 * independent of contact-force laws and structural time integration.
 */
class NarrowPhase {
public:
  /**
   * @brief Destroy the interface polymorphically.
   */
  virtual ~NarrowPhase() = 0;
};

inline NarrowPhase::~NarrowPhase() = default;

} // namespace ELFF::Contact
