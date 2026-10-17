#include <gtest/gtest.h>

#include "elff/c/models/ibm/IBRuntime.h"
#include "elff/c/models/ibm/IBVelocityCoupled.h"
#include "elff/c/models/spring/spring.h"

namespace {

TEST(IBSpringCAPITest, RuntimeOwnsCircleAndEllipse)
{
  ib_runtime_t runtime = ib_runtime_new();
  ASSERT_NE(runtime, nullptr);
  ib_spring_circle_t circle =
    ib_spring_circle_create_circle(runtime, 16, 1.0, 0.5, 1.0, 2.0);
  ib_spring_circle_t ellipse =
    ib_spring_circle_create_ellipse(runtime, 24, 2.0, 0.5, 0.25, -1.0, 3.0);
  ASSERT_NE(circle, nullptr);
  ASSERT_NE(ellipse, nullptr);
  EXPECT_EQ(ib_velocity_coupled_get_number_of_nodes(circle), 16);
  EXPECT_EQ(ib_velocity_coupled_get_number_of_nodes(ellipse), 24);

  ib_model_view_array_t views = {};
  ASSERT_EQ(ib_runtime_get_model_views(runtime, &views), 0);
  ASSERT_EQ(views.count, 2U);
  EXPECT_EQ(views.data[0].handle, circle);
  EXPECT_EQ(views.data[1].handle, ellipse);
  EXPECT_EQ(views.data[0].kind, IB_MODEL_VIEW_VELOCITY_COUPLED);
  ib_runtime_delete(runtime);
}

} // namespace
