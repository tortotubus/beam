#include "elff/contact/NarrowPhaseSegmentDistance.hpp"
#include "elff/contact/BroadPhaseLinkedCell.hpp"
#include "elff/general/error.hpp"
#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <random>

using namespace ELFF::Contact;
namespace {
using ELFF::real_t;
const NarrowPhaseSegmentDistance narrow;
#ifdef ELFF_USE_EXCEPTIONS
#define EXPECT_SEGMENT_ERROR(statement) EXPECT_THROW(statement, ELFF::ErrorException)
#else
#define EXPECT_SEGMENT_ERROR(statement) EXPECT_DEATH({ statement; }, "BEAM abort:")
#endif

TEST(NarrowPhaseSegmentDistanceTest, InteriorCrossingAndSkew)
{
  const LineSegment a{{-1,0,0},{1,0,0}};
  for (real_t z : {0., 2.}) {
    const LineSegment b{{0,-1,z},{0,1,z}};
    const auto result = narrow.closest_points(a,b);
    EXPECT_DOUBLE_EQ(result.distance_squared, z*z);
    EXPECT_DOUBLE_EQ(result.parameter_a, 0.5);
    EXPECT_DOUBLE_EQ(result.parameter_b, 0.5);
    EXPECT_EQ(result.point_a, (std::array<real_t,3>{0,0,0}));
    EXPECT_EQ(result.point_b, (std::array<real_t,3>{0,0,z}));
    EXPECT_TRUE(narrow.within_clearance(a,b,z));
    if (z > 0) EXPECT_FALSE(narrow.within_clearance(a,b,std::nextafter(z,0.)));
  }
}

TEST(NarrowPhaseSegmentDistanceTest, EndpointAndParallelMinima)
{
  const LineSegment a{{0,0,0},{1,0,0}};
  const LineSegment endpoint{{2,1,0},{3,2,0}};
  auto result = narrow.closest_points(a,endpoint);
  EXPECT_DOUBLE_EQ(result.distance_squared, 2);
  EXPECT_DOUBLE_EQ(result.parameter_a, 1);
  EXPECT_DOUBLE_EQ(result.parameter_b, 0);
  EXPECT_DOUBLE_EQ(narrow.closest_points(a,{{0.5,2,0},{2,2,0}}).distance_squared,4);
  EXPECT_DOUBLE_EQ(narrow.closest_points(a,{{2,2,0},{0.5,2,0}}).distance_squared,4);
  EXPECT_DOUBLE_EQ(narrow.closest_points(a,{{0.5,0,0},{2,0,0}}).distance_squared,0);
  EXPECT_DOUBLE_EQ(narrow.closest_points(a,{{2,0,0},{3,0,0}}).distance_squared,1);
  EXPECT_TRUE(narrow.within_clearance(a,{{1,0,0},{2,0,0}},0));
}

TEST(NarrowPhaseSegmentDistanceTest, DegenerateAndNearlyParallelSegments)
{
  const LineSegment point{{0.5,1,0},{0.5,1,0}}, line{{0,0,0},{1,0,0}};
  EXPECT_DOUBLE_EQ(narrow.closest_points(point,line).distance_squared,1);
  EXPECT_DOUBLE_EQ(narrow.closest_points(line,point).distance_squared,1);
  EXPECT_DOUBLE_EQ(narrow.closest_points(point,point).distance_squared,0);
  EXPECT_DOUBLE_EQ(narrow.closest_points(point,{{0.5,3,0},{0.5,3,0}}).distance_squared,4);
  // A threshold-based parallel shortcut would incorrectly choose an endpoint.
  const auto result = narrow.closest_points(line,{{0,-1e-12,0},{1,1e-12,0}});
  EXPECT_NEAR(result.parameter_a,0.5,1e-14);
  EXPECT_NEAR(result.parameter_b,0.5,1e-14);
  EXPECT_NEAR(result.distance_squared,0,1e-40);
}

TEST(NarrowPhaseSegmentDistanceTest, SymmetryWitnessesAndOptimality)
{
  std::mt19937 generator(42);
  std::uniform_real_distribution<real_t> coordinate(-5,5);
  for (int trial = 0; trial < 300; ++trial) {
    LineSegment a,b;
    for (auto* point : {&a.start,&a.end,&b.start,&b.end})
      for (auto& value : *point) value = coordinate(generator);
    const auto result = narrow.closest_points(a,b);
    EXPECT_NEAR(result.distance_squared,narrow.closest_points(b,a).distance_squared,1e-12);
    EXPECT_NEAR(result.distance_squared,narrow.closest_points({a.end,a.start},b).distance_squared,1e-12);
    real_t squared = 0, ds = 0, dt = 0;
    EXPECT_GE(result.parameter_a,0); EXPECT_LE(result.parameter_a,1);
    EXPECT_GE(result.parameter_b,0); EXPECT_LE(result.parameter_b,1);
    for (int d = 0; d < 3; ++d) {
      EXPECT_NEAR(result.point_a[d],a.start[d]+result.parameter_a*(a.end[d]-a.start[d]),1e-14);
      EXPECT_NEAR(result.point_b[d],b.start[d]+result.parameter_b*(b.end[d]-b.start[d]),1e-14);
      const real_t delta = result.point_a[d]-result.point_b[d];
      squared += delta*delta;
      ds += delta*(a.end[d]-a.start[d]);
      dt -= delta*(b.end[d]-b.start[d]);
    }
    EXPECT_NEAR(squared,result.distance_squared,1e-12);
    // First-order conditions are sufficient for this convex box-constrained problem.
    for (const auto entry : {std::pair{result.parameter_a,ds},std::pair{result.parameter_b,dt}}) {
      if (entry.first == 0) EXPECT_GE(entry.second,-1e-12);
      else if (entry.first == 1) EXPECT_LE(entry.second,1e-12);
      else EXPECT_NEAR(entry.second,0,1e-12);
    }
  }
}

TEST(NarrowPhaseSegmentDistanceTest, GridCandidatesSupportClearanceRejection)
{
  const LineSegment a{{1,1,1},{3,1,1}}, b{{2,1.25,1},{2,2,1}};
  BroadPhaseLinkedCell grid({0,0,0},1,{5,5,5});
  // Expand one AABB by the full centerline clearance.
  grid.insert(0,{{0.5,0.5,0.5},{3.5,1.5,1.5}});
  grid.insert(1,{{2,1.25,1},{2,2,1}});
  EXPECT_EQ(grid.candidate_pairs(),(std::vector<ShapePair>{{0,1}}));
  EXPECT_TRUE(narrow.within_clearance(a,b,0.5));
  EXPECT_FALSE(narrow.within_clearance(a,b,0.125));
}

TEST(NarrowPhaseSegmentDistanceTest, SmallAndLargeCoordinateScales)
{
  for (real_t scale : {1e-150, 1e150}) {
    const LineSegment a{{-scale,0,0},{scale,0,0}};
    const LineSegment b{{0,-scale,scale},{0,scale,scale}};
    const auto result = narrow.closest_points(a,b);
    EXPECT_NEAR(result.distance_squared/(scale*scale),1,1e-14);
    EXPECT_DOUBLE_EQ(result.parameter_a,0.5);
    EXPECT_DOUBLE_EQ(result.parameter_b,0.5);
    EXPECT_FALSE(narrow.within_clearance(a,b,scale/2));
    EXPECT_TRUE(narrow.within_clearance(a,b,2*scale));
  }
}

TEST(NarrowPhaseSegmentDistanceTest, RejectsInvalidInputs)
{
  const LineSegment a{{0,0,0},{1,0,0}};
  const auto nan = std::numeric_limits<real_t>::quiet_NaN();
  EXPECT_SEGMENT_ERROR(narrow.closest_points(a,{{nan,0,0},{1,0,0}}));
  EXPECT_SEGMENT_ERROR(narrow.within_clearance(a,a,-1));
  EXPECT_SEGMENT_ERROR(narrow.within_clearance(a,a,nan));
}
#undef EXPECT_SEGMENT_ERROR
} // namespace
