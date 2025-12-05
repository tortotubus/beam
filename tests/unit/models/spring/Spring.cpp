#include <gtest/gtest.h>

#include "models/spring/Spring.hpp"

using namespace beam;

// TEST(IBSpringCircleTest, ConstructionAndPoints)
// {
//     int N = 8;
//     real_t K = 1.0;
//     real_t radius = 2.0;
//     real_t cx = 0.0, cy = 0.0;
//     IBSpringCircle circle(N, K, radius, cx, cy);

//     // Check number of points
//     EXPECT_EQ(circle.GetNumberOfPoints(), N);

//     // Check that all points are on the circle (distance from center == radius)
//     auto& mesh = circle.GetCurrent();
//     auto& points = mesh.GetPoints();
//     for (int i = 0; i < N; ++i) {
//         real_t x = points[i](0);
//         real_t y = points[i](1);
//         real_t dist = std::sqrt((x - cx)*(x - cx) + (y - cy)*(y - cy));
//         EXPECT_NEAR(dist, radius, 1e-6);
//     }
// }

// TEST(IBSpringCircleTest, GetMidpointEllipseZeroVelocityProducesNonZeroForces)
// {
//     int N = 12;
//     real_t K = 1.0;
//     real_t rx = 2.0, ry = 1.0;
//     real_t cx = 0.0, cy = 0.0;
//     IBSpringCircle ellipse(N, K, rx, ry, cx, cy);

//     // Zero velocity array
//     Array<fem::Vertex> velocity(N);
//     for (int i = 0; i < N; ++i) {
//         velocity[i](0) = 0.0;
//         velocity[i](1) = 0.0;
//         // If Vertex is 3D, set z to 0 as well (assume up to 3D)
//         try { velocity[i](2) = 0.0; } catch (...) {}
//     }
//     real_t dt = 0.1;

//     // GetMidpoint should not move points, but should compute forces
//     IBStructureMesh& mesh = ellipse.GetMidpoint(velocity, dt);
//     auto& forces = mesh.GetForces();

//     // At least one force should be nonzero (ellipse is not equilibrium)
//     bool any_nonzero = false;
//     for (int i = 0; i < N; ++i) {
//         real_t fx = forces[i](0);
//         real_t fy = forces[i](1);
//         if (std::abs(fx) > 1e-8 || std::abs(fy) > 1e-8) {
//             any_nonzero = true;
//             break;
//         }
//     }
//     EXPECT_TRUE(any_nonzero) << "Forces should be nonzero for ellipse at rest";
// }