#pragma once

#include "elff/contact/BroadPhase.hpp"
#include <algorithm>
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

  /** @brief Apply eligibility before emitting, sorting or deduplicating pairs.
   * The predicate receives inserted shape IDs (in arbitrary order), can be
   * evaluated repeatedly for the same pair, and must not modify the grid.
   */
  template <class Predicate> std::vector<ShapePair> candidate_pairs_if(Predicate eligible) const {
    std::vector<ShapePair> pairs;
    for (const auto& cell : cells_)
      for (std::size_t i = 0; i < cell.shapes.size(); ++i)
        for (std::size_t j = i + 1; j < cell.shapes.size(); ++j) {
          const auto a = cell.shapes[i], b = cell.shapes[j];
          if (a != b && eligible(a, b))
            pairs.emplace_back(std::minmax(a, b));
        }
    std::sort(pairs.begin(), pairs.end());
    pairs.erase(std::unique(pairs.begin(), pairs.end()), pairs.end());
    return pairs;
  }

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
  struct Cell {
    std::vector<ShapeId> shapes;
  };

  std::array<real_t, 3> origin_;          /**< @brief Lower domain corner. */
  std::array<real_t, 3> upper_;           /**< @brief Closed upper domain corner. */
  real_t cell_width_;                   /**< @brief Uniform cell width. */
  std::array<std::size_t, 3> dimensions_; /**< @brief Cell counts in each direction. */
  std::vector<Cell> cells_; /**< @brief Dense storage, indexed by x + nx * (y + ny * z). */
};

} // namespace ELFF::Contact
