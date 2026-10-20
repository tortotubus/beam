#pragma once

namespace ELFF::Contact {

/**
 * @brief Orchestrates contact detection, force evaluation and load assembly.
 *
 * Implementations should coordinate geometry adapters, candidate generation,
 * pair filtering, narrow-phase queries, contact laws and load mapping. The
 * system should manage trial and committed contact history while leaving
 * structural state advancement to the owning solver.
 */
class ContactSystem {
public:
  /**
   * @brief Destroy the interface polymorphically.
   */
  virtual ~ContactSystem() = 0;
};

inline ContactSystem::~ContactSystem() = default;

} // namespace ELFF::Contact
