#include "elff/contact/BroadPhaseLinkedCell.hpp"
#include "elff/general/error.hpp"
#include <algorithm>
#include <cmath>
#include <limits>

namespace ELFF::Contact {

BroadPhaseLinkedCell::BroadPhaseLinkedCell(std::array<real_t, 3> origin,
    real_t cell_width, std::array<std::size_t, 3> dimensions)
  : origin_(origin), cell_width_(cell_width), dimensions_(dimensions)
{
  if (!std::isfinite(cell_width_) || cell_width_ <= 0)
    ELFF_ABORT("Linked-cell width must be finite and positive.");
  std::size_t count = 1;
  for (std::size_t d = 0; d < 3; ++d) {
    if (!dimensions_[d] || dimensions_[d] > cells_.max_size() / count)
      ELFF_ABORT("Invalid or overflowing linked-cell dimensions.");
    count *= dimensions_[d];
    upper_[d] = origin_[d] + cell_width_ * dimensions_[d];
    if (!std::isfinite(origin_[d]) || !std::isfinite(upper_[d]) ||
        upper_[d] <= origin_[d])
      ELFF_ABORT("Linked-cell domain must be finite with positive extent.");
  }
  cells_.resize(count);
}

void BroadPhaseLinkedCell::clear()
{
  for (auto& cell : cells_) cell.shapes.clear();
}

void BroadPhaseLinkedCell::insert(ShapeId id, const AABB& bounds)
{
  std::array<std::size_t, 3> first, last;
  for (std::size_t d = 0; d < 3; ++d) {
    if (!std::isfinite(bounds.min[d]) || !std::isfinite(bounds.max[d]) ||
        bounds.min[d] > bounds.max[d] || bounds.min[d] < origin_[d] ||
        bounds.max[d] > upper_[d])
      ELFF_ABORT("Linked-cell bounds must be finite, ordered and inside the domain.");
    auto index = [&](real_t x, bool lower) {
      // Round outward conservatively at cell faces, then clip the closed
      // domain boundary to its adjacent cell. Use wider arithmetic for offsets.
      long double q = (static_cast<long double>(x) - origin_[d]) / cell_width_;
      q = std::floor(std::nextafter(q, lower
          ? -std::numeric_limits<long double>::infinity()
          : std::numeric_limits<long double>::infinity()));
      if (q <= 0) return std::size_t{0};
      if (q >= static_cast<long double>(dimensions_[d] - 1))
        return dimensions_[d] - 1;
      return static_cast<std::size_t>(q);
    };
    first[d] = index(bounds.min[d], true);
    last[d] = index(bounds.max[d], false);
  }
  // Validate every coordinate before modifying any cells.
  for (std::size_t z = first[2]; z <= last[2]; ++z)
    for (std::size_t y = first[1]; y <= last[1]; ++y)
      for (std::size_t x = first[0]; x <= last[0]; ++x)
        cells_[x + dimensions_[0] * (y + dimensions_[1] * z)].shapes.push_back(id);
}

std::vector<ShapePair> BroadPhaseLinkedCell::candidate_pairs() const
{
  std::vector<ShapePair> pairs;
  for (const auto& cell : cells_)
    for (std::size_t i = 0; i < cell.shapes.size(); ++i)
      for (std::size_t j = i + 1; j < cell.shapes.size(); ++j)
        if (cell.shapes[i] != cell.shapes[j])
          pairs.emplace_back(std::minmax(cell.shapes[i], cell.shapes[j]));
  std::sort(pairs.begin(), pairs.end());
  pairs.erase(std::unique(pairs.begin(), pairs.end()), pairs.end());
  return pairs;
}

std::vector<ShapePair>
BroadPhaseLinkedCell::candidate_pairs(ShapeId first_candidate) const
{
  std::vector<ShapePair> pairs;
  for (const auto& cell : cells_)
    for (const auto candidate : cell.shapes) {
      if (candidate < first_candidate) continue;
      for (const auto obstacle : cell.shapes)
        if (obstacle < first_candidate)
          pairs.emplace_back(obstacle, candidate);
    }
  std::sort(pairs.begin(), pairs.end());
  pairs.erase(std::unique(pairs.begin(), pairs.end()), pairs.end());
  return pairs;
}

} // namespace ELFF::Contact
