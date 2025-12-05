#include <gtest/gtest.h>

extern "C" {

#include "basilisk/interface/ibm/immersedboundary.h"
#include "basilisk/interface/ibm/spring.h"
}

#include <vector>
#include <cmath>


TEST(IBSpringCircleCInterfaceTest, GetMidpointEllipseZeroVelocityProducesNonZeroForces)
{
    // Create the spring object as an ellipse using the C interface
    int N = 12;
    double K = 1.0;
    double rx = 2.0, ry = 1.0;
    double cx = 0.0, cy = 0.0;
    ib_velocity_structure_model_t handle = ib_spring_circle_create_ellipse(N, K, rx, ry, cx, cy);
    ASSERT_NE(handle, nullptr);

    // Zero velocity array
    std::vector<vertex_t> velocity(N);
    for (int i = 0; i < N; ++i) {
        velocity[i].x = 0.0;
        velocity[i].y = 0.0;
        velocity[i].z = 0.0;
    }
    double dt = 0.1;

    // Call the C interface for midpoint
    ib_structure_mesh_t mesh = ib_velocity_structure_model_get_midpoint(handle, velocity.data(), N, dt);

    // Check that at least one force is nonzero
    bool any_nonzero = false;
    for (int i = 0; i < mesh.n; ++i) {
        if (std::abs(mesh.forces[i].x) > 1e-8 || std::abs(mesh.forces[i].y) > 1e-8) {
            any_nonzero = true;
            break;
        }
    }
    EXPECT_TRUE(any_nonzero) << "Forces should be nonzero for ellipse at rest (C interface)";

    // Free memory allocated by the C interface
    ib_structure_mesh_free(&mesh);
    ib_spring_circle_destroy(handle);
}
