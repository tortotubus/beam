#include "spring.h"

#include "models/spring/Spring.hpp"

using namespace beam;

extern "C"
{

  ib_spring_circle_t ib_spring_circle_create_circle(int npoints,
                                             double K,
                                             double radius,
                                             double center_x,
                                             double center_y)
  {
    return new IBSpringCircle(npoints, K, radius, center_x, center_y);
  }

  ib_spring_circle_t ib_spring_circle_create_ellipse(int npoints,
                                             double K,
                                             double radius_x,
                                             double radius_y,
                                             double center_x,
                                             double center_y)
  {
    return new IBSpringCircle(
      npoints, K, radius_x, radius_y, center_x, center_y);
  }

  void ib_spring_circle_destroy(ib_spring_circle_t handle)
  {
    delete reinterpret_cast<IBSpringCircle*>(handle);
  }
}