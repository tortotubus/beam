#include "elff/contact/BroadPhaseLinkedCell.hpp"
#include "elff/general/error.hpp"
#include <gtest/gtest.h>
#include <algorithm>
#include <limits>
#include <random>

using namespace ELFF::Contact;
namespace {
#ifdef ELFF_USE_EXCEPTIONS
#define EXPECT_GRID_ERROR(statement) EXPECT_THROW(statement, ELFF::ErrorException)
#else
#define EXPECT_GRID_ERROR(statement) EXPECT_DEATH({ statement; }, "BEAM abort:")
#endif

TEST(BroadPhaseLinkedCellTest, DenseCellsDeduplicateAndClear)
{
  BroadPhaseLinkedCell grid({-2, -2, -2}, 1, {4, 4, 4});
  BroadPhase& broad = grid;
  broad.insert(9, {{-1.8, -1.8, -1.8}, {0.8, 0.8, 0.8}});
  broad.insert(2, {{-1.5, -1.5, -1.5}, {0.5, 0.5, 0.5}});
  broad.insert(2, {{-1.5, -1.5, -1.5}, {0.5, 0.5, 0.5}});
  broad.insert(7, {{1.2, 1.2, 1.2}, {1.8, 1.8, 1.8}});
  EXPECT_EQ(broad.candidate_pairs(), (std::vector<ShapePair>{{2, 9}}));
  broad.clear();
  EXPECT_TRUE(broad.candidate_pairs().empty());
  broad.insert(7, {{2, 2, 2}, {2, 2, 2}});
  broad.insert(3, {{1.9, 1.9, 1.9}, {2, 2, 2}});
  EXPECT_EQ(broad.candidate_pairs(), (std::vector<ShapePair>{{3, 7}}));
}

TEST(BroadPhaseLinkedCellTest, TouchingFacesEdgesAndCornersAreCandidates)
{
  BroadPhaseLinkedCell grid({0, 0, 0}, 1, {3, 3, 3});
  grid.insert(0, {{0, 0, 0}, {1, 1, 1}});
  grid.insert(1, {{1, 0, 0}, {2, 1, 1}});
  grid.insert(2, {{1, 1, 0}, {2, 2, 1}});
  grid.insert(3, {{1, 1, 1}, {2, 2, 2}});
  grid.insert(4, {{0, 0, 0}, {0, 0, 0}});
  const auto pairs = grid.candidate_pairs();
  for (ShapeId id = 1; id <= 4; ++id)
    EXPECT_TRUE(std::binary_search(pairs.begin(), pairs.end(), ShapePair{0, id}));
}

TEST(BroadPhaseLinkedCellTest, ContainsEveryAllPairsOverlap)
{
  BroadPhaseLinkedCell grid({-3, -3, -3}, 0.3, {20, 20, 20});
  std::mt19937 random(1234);
  std::uniform_real_distribution<double> position(-3, 2);
  std::uniform_real_distribution<double> extent(0, 1);
  std::vector<AABB> bounds(150);
  for (ShapeId id = 0; id < bounds.size(); ++id) {
    for (int d = 0; d < 3; ++d) {
      bounds[id].min[d] = position(random);
      bounds[id].max[d] = bounds[id].min[d] + extent(random);
    }
    grid.insert(id, bounds[id]);
  }
  const auto pairs = grid.candidate_pairs();
  EXPECT_TRUE(std::is_sorted(pairs.begin(), pairs.end()));
  EXPECT_EQ(std::adjacent_find(pairs.begin(), pairs.end()), pairs.end());
  for (ShapeId a = 0; a < bounds.size(); ++a)
    for (ShapeId b = a + 1; b < bounds.size(); ++b) {
      bool overlap = true;
      for (int d = 0; d < 3; ++d)
        overlap &= bounds[a].min[d] <= bounds[b].max[d] &&
                   bounds[b].min[d] <= bounds[a].max[d];
      if (overlap)
        EXPECT_TRUE(std::binary_search(pairs.begin(), pairs.end(), ShapePair{a, b}));
    }
}

TEST(BroadPhaseLinkedCellTest, GroupPairsMatchFilteredAllPairs)
{
  BroadPhaseLinkedCell grid({0, 0, 0}, 1, {4, 4, 4});
  std::mt19937 random(42);
  std::uniform_real_distribution<double> position(0, 3);
  // Scrambled insertion order and repeated IDs must not affect group membership.
  for (ShapeId n = 0; n < 80; ++n) {
    const ShapeId id = (n * 17) % 40;
    const std::array<ELFF::real_t, 3> lo{
        position(random), position(random), position(random)};
    grid.insert(id, {lo, {lo[0] + 1, lo[1] + 1, lo[2] + 1}});
  }
  grid.insert(0, {{0, 0, 0}, {1, 1, 1}});
  grid.insert(40, {{1, 1, 1}, {1, 1, 1}});
  grid.insert(41, {{4, 4, 4}, {4, 4, 4}});
  grid.insert(1, {{3, 3, 3}, {4, 4, 4}});
  const auto all = grid.candidate_pairs();
  for (ShapeId split : {0, 1, 20, 40, 42}) {
    std::vector<ShapePair> expected;
    for (const auto& pair : all)
      if (pair.first < split && pair.second >= split)
        expected.push_back(pair);
    EXPECT_EQ(grid.candidate_pairs(split), expected);
  }
  grid.clear();
  EXPECT_TRUE(grid.candidate_pairs(20).empty());
}

TEST(BroadPhaseLinkedCellTest, RejectsInvalidGridAndBounds)
{
  const auto nan = std::numeric_limits<ELFF::real_t>::quiet_NaN();
  EXPECT_GRID_ERROR(BroadPhaseLinkedCell({0, 0, 0}, 0, {1, 1, 1}));
  EXPECT_GRID_ERROR(BroadPhaseLinkedCell({0, 0, 0}, 1, {0, 1, 1}));
  EXPECT_GRID_ERROR(BroadPhaseLinkedCell({nan, 0, 0}, 1, {1, 1, 1}));
  EXPECT_GRID_ERROR(BroadPhaseLinkedCell({0, 0, 0}, 1,
      {std::numeric_limits<std::size_t>::max(), 2, 2}));
  BroadPhaseLinkedCell grid({0, 0, 0}, 1, {2, 2, 2});
  grid.insert(0, {{0, 0, 0}, {2, 2, 2}});
  EXPECT_GRID_ERROR(grid.insert(1, {{0, 0, -0.1}, {1, 1, 1}}));
  EXPECT_GRID_ERROR(grid.insert(1, {{0, 0, 0}, {1, 1, 2.1}}));
  EXPECT_GRID_ERROR(grid.insert(1, {{1, 0, 0}, {0, 1, 1}}));
  EXPECT_GRID_ERROR(grid.insert(1, {{0, 0, 0}, {1, 1, nan}}));
  EXPECT_TRUE(grid.candidate_pairs().empty());
}
#undef EXPECT_GRID_ERROR
} // namespace
