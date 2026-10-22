#pragma once

#include "elff/contact/BroadPhase.hpp"
#include <cstddef>

namespace ELFF::Contact {

/**
 * @brief Dense uniform grid for conservative broad-phase contact detection.
 *
 * Each shape is inserted into every cell touched by its closed AABB. Candidate
 * pairs are generated within cells and deduplicated; no neighboring-cell
 * traversal is required. The caller controls clear() and insertion explicitly.
 * Bounds outside the fixed domain are rejected through ELFF_ABORT.
 */
class BroadPhaseLinkedCell : public BroadPhase {
public:
  /**
   * @brief Allocate a grid with equal cell width in all three directions.
   * @param origin Lower corner of the domain.
   * @param cell_width Positive finite cell width.
   * @param dimensions Positive cell counts in x, y and z.
   * The upper domain corner is origin plus cell_width times dimensions.
   */
  BroadPhaseLinkedCell(std::array<real_t, 3> origin, real_t cell_width,
                      std::array<std::size_t, 3> dimensions);

  /** @brief Empty all cells while retaining their allocated storage. */
  void clear() override;

  /**
   * @brief Insert finite, ordered bounds fully contained in the closed domain.
   * @param id Stable shape identifier. Repeated IDs never form self-pairs.
   * @param bounds World-space bounds including radius and contact-search margin.
   * Zero-width bounds and bounds touching the domain boundary are accepted.
   */
  void insert(ShapeId id, const AABB& bounds) override;

  /** @brief Return sorted unique candidate pairs; geometric false positives are allowed. */
  std::vector<ShapePair> candidate_pairs() const override;

  /**
   * @brief Return sorted unique pairs across two disjoint groups of shape IDs.
   * @param first_candidate IDs below this value are obstacles; all others are candidates.
   *
   * Only obstacle-candidate pairs are generated. Pairs within either group are
   * never materialized, avoiding unnecessary work during placement queries.
   */
  std::vector<ShapePair> candidate_pairs(ShapeId first_candidate) const;

private:
  /** @brief Shape identifiers occupying a single grid cell. */
  struct Cell { std::vector<ShapeId> shapes; };

  std::array<real_t, 3> origin_; /**< @brief Lower domain corner. */
  std::array<real_t, 3> upper_; /**< @brief Closed upper domain corner. */
  real_t cell_width_; /**< @brief Uniform cell width. */
  std::array<std::size_t, 3> dimensions_; /**< @brief Cell counts in each direction. */
  std::vector<Cell> cells_; /**< @brief Dense storage, indexed by x + nx * (y + ny * z). */
};

} // namespace ELFF::Contact
