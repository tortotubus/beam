#pragma once

#include "elff/config/config.hpp"
#include <array>
#include <cstdint>
#include <utility>
#include <vector>

namespace ELFF::Contact {

using ShapeId = std::uint64_t;
using ShapePair = std::pair<ShapeId, ShapeId>;

/** @brief Closed world-space bounds, including any physical radius and search margin. */
struct AABB {
  std::array<real_t, 3> min;
  std::array<real_t, 3> max;
};

/**
 * @brief Identifies candidate pairs of collision shapes.
 *
 * Implementations should use spatial indexing or bounding volumes to return
 * potentially interacting shape pairs for NarrowPhase. Candidate generation
 * may include false positives, but must retain every pair within the contact
 * search range.
 */
class BroadPhase {
public:
  /**
   * @brief Destroy the interface polymorphically.
   */
  virtual ~BroadPhase() = 0;

  /** @brief Remove all inserted shapes while retaining the spatial configuration. */
  virtual void clear() = 0;

  /** @brief Insert a shape's conservative world-space bounds. */
  virtual void insert(ShapeId id, const AABB& bounds) = 0;

  /** @brief Return unique candidate pairs, with the smaller ID first, in sorted order. */
  virtual std::vector<ShapePair> candidate_pairs() const = 0;
};

inline BroadPhase::~BroadPhase() = default;

} // namespace ELFF::Contact
