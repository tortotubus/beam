#pragma once

namespace ELFF::Contact {

/**
 * @brief Determines which candidate shape pairs may interact.
 *
 * Implementations should apply model-specific interaction rules before
 * narrow-phase detection, such as excluding adjacent rod links or disabled
 * material pairs. Filtering determines whether a pair is eligible for contact;
 * it does not determine whether the shapes geometrically intersect.
 */
class CollisionFilter {
public:
  /**
   * @brief Destroy the interface polymorphically.
   */
  virtual ~CollisionFilter() = 0;
};

inline CollisionFilter::~CollisionFilter() = default;

} // namespace ELFF::Contact
