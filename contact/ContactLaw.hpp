#pragma once

namespace ELFF::Contact {

/**
 * @brief Evaluates contact forces and moments from geometry and motion.
 *
 * Implementations should evaluate normal, tangential and optional rolling
 * responses from contact geometry, relative motion, material parameters and
 * committed ContactHistory. Results should provide forces, moments and trial
 * history updates independently of the collision search and structural solver.
 */
class ContactLaw {
public:
  /**
   * @brief Destroy the interface polymorphically.
   */
  virtual ~ContactLaw() = 0;
};

inline ContactLaw::~ContactLaw() = default;

} // namespace ELFF::Contact
