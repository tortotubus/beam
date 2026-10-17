#include "elff/c/models/spring/spring.h"

#include "elff/models/spring/IBSpringCircle.hpp"
#include "elff/models/ibm/IBRuntime.hpp"

#include <cmath>
#include <memory>

using namespace ELFF;
using namespace ELFF::Models;

extern "C"
{

  ib_spring_circle_t ib_spring_circle_create_circle(ib_runtime_t runtime,
                                             int npoints,
                                             double K,
                                             double radius,
                                             double center_x,
                                             double center_y)
  {
    if (!runtime || npoints < 3 || !std::isfinite(K) ||
        !std::isfinite(radius) || radius <= 0.0)
      return nullptr;
    try {
      auto spring = std::make_unique<IBSpringCircle>(
        npoints, K, radius, center_x, center_y);
      return reinterpret_cast<Models::IBRuntime*>(runtime)->adopt_owned_model(
        Models::IBRuntimeModelType::SpringCircle, std::move(spring));
    } catch (...) {
      return nullptr;
    }
  }

  ib_spring_circle_t ib_spring_circle_create_ellipse(ib_runtime_t runtime,
                                             int npoints,
                                             double K,
                                             double radius_x,
                                             double radius_y,
                                             double center_x,
                                             double center_y)
  {
    if (!runtime || npoints < 3 || !std::isfinite(K) ||
        !std::isfinite(radius_x) || !std::isfinite(radius_y) ||
        radius_x <= 0.0 || radius_y <= 0.0)
      return nullptr;
    try {
      auto spring = std::make_unique<IBSpringCircle>(
        npoints, K, radius_x, radius_y, center_x, center_y);
      return reinterpret_cast<Models::IBRuntime*>(runtime)->adopt_owned_model(
        Models::IBRuntimeModelType::SpringCircle, std::move(spring));
    } catch (...) {
      return nullptr;
    }
  }
}
