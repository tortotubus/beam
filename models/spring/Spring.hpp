#pragma once

#include "../immersedboundary.hpp"

#include <cmath>

namespace beam {
class IBSpringCircle : public IBStructureModel
{
protected:
  real_t K;
  real_t radius_x, radius_y;
  real_t dtheta;

public:
  IBSpringCircle(int NumberOfPoints,
                 real_t K,
                 real_t radius,
                 real_t center_x,
                 real_t center_y)
    : IBStructureModel(2, NumberOfPoints)
    , K(K)
    , radius_x(radius)
    , radius_y(radius)
    , dtheta((2 * M_1_PI) / NumberOfPoints)
  {
    Array<fem::Vertex> circle_points =
      circle(NumberOfPoints, radius, center_x, center_y);

    for (int i = 0; i < NumberOfPoints; ++i) {
      mesh.GetPoints()[i] = circle_points[i];
    }
  };

  IBSpringCircle(int NumberOfPoints,
                 real_t K,
                 real_t radius_x,
                 real_t radius_y,
                 real_t center_x,
                 real_t center_y)
    : IBStructureModel(2, NumberOfPoints)
    , K(K)
    , radius_x(radius_x)
    , radius_y(radius_y)
    , dtheta((2 * M_1_PI) / NumberOfPoints)
  {
    Array<fem::Vertex> ellipse_point =
      ellipse(NumberOfPoints, radius_x, radius_y, center_x, center_y);

    for (int i = 0; i < NumberOfPoints; ++i) {
      mesh.GetPoints()[i] = ellipse_point[i];
    }
  };

protected:
  void ComputeMidpointForces() override
  {
    Array<fem::Vertex>& forces = mesh_next.GetForces();
    Array<fem::Vertex>& midpoints = mesh_next.GetPoints();

    int N = GetNumberOfPoints();
    real_t inv_dtheta2 = 1. / (dtheta * dtheta);

    for (int i = 0; i < N; i++) {
      int ip = (i - 1 + N) % N;
      int in = (i + 1) % N;
      forces[i](0) = K * (midpoints[ip](0) + midpoints[in](0) - 2 * midpoints[i](0)) * inv_dtheta2;
      forces[i](1) = K * (midpoints[ip](1) + midpoints[in](1) - 2 * midpoints[i](1)) * inv_dtheta2;
    }
  }

private:
  static Array<fem::Vertex> circle(int N,
                                   real_t radius,
                                   real_t center_x,
                                   real_t center_y)
  {
    Array<fem::Vertex> points(N);
    real_t dtheta = (2. * M_PI) / N;
    for (int pi = 0; pi < N; pi++) {
      points[pi](0) = radius * cos(dtheta * pi) + center_x;
      points[pi](1) = radius * sin(dtheta * pi) + center_y;
    }
    return points;
  }

  static Array<fem::Vertex> ellipse(int N,
                                    real_t radius_x,
                                    real_t radius_y,
                                    real_t center_x,
                                    real_t center_y)
  {
    Array<fem::Vertex> points(N);
    real_t dtheta = (2. * M_PI) / N;
    for (int pi = 0; pi < N; pi++) {
      points[pi](0) = radius_x * cos(dtheta * pi) + center_x;
      points[pi](1) = radius_y * sin(dtheta * pi) + center_y;
    }
    return points;
  }
};
} // namespace beam