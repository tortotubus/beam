#pragma once

namespace ELFF::Contact {

/**
 * @brief Manages persistent and trial constitutive contact history.
 *
 * Implementations should associate friction displacement, previous normals
 * and other constitutive memory with persistent contact identities. Trial
 * evaluations must use committed history without accumulating changes across
 * nonlinear iterations; accepted steps commit updates and rejected steps
 * discard them.
 */
class ContactHistory {
public:
  /**
   * @brief Destroy the interface polymorphically.
   */
  virtual ~ContactHistory() = 0;
};

inline ContactHistory::~ContactHistory() = default;

} // namespace ELFF::Contact
