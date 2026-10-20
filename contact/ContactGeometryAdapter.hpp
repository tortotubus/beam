#pragma once

namespace ELFF::Contact {

/**
 * @brief Exposes collision geometry and kinematics from a model state.
 *
 * Implementations should expose collision shapes, their positions and
 * orientations, and contact-point velocities from a model's current or trial
 * state. The structural solver retains ownership of that state and its time
 * integration.
 */
class ContactGeometryAdapter {
public:
  /**
   * @brief Destroy the interface polymorphically.
   */
  virtual ~ContactGeometryAdapter() = 0;
};

inline ContactGeometryAdapter::~ContactGeometryAdapter() = default;

} // namespace ELFF::Contact
