#pragma once

#include "../immersedboundary.hpp"

namespace beam {
class IBSpringCircle : public IBStructureModel
{
protected:
  real_t K;
  real_t radius;
  real_t dtheta;
public:
  IBSpringCircle(int NumberOfPoints,
                     real_t K,
                     real_t radius,
                     real_t center_x,
                     real_t center_y)
    : IBStructureModel(2, NumberOfPoints)
    , K(K)
    , radius(radius)
    , dtheta((2*M_1_PI) / NumberOfPoints)
  {
    Array<real_t> circle_points =
      circle(NumberOfPoints, radius, center_x, center_y);
    mesh.GetPoints().Copy(circle_points);
  };
protected:
  void ComputeMidpointForces() override {
    Array<real_t> &forces = mesh_next.GetForces();
    Array<real_t> &midpoints = mesh_next.GetPoints();

    int N = GetNumberOfPoints();
    real_t inv_dtheta2 = 1./(dtheta*dtheta);

    for (int i = 0; i < N; i++) {
      int ip = (i-1) % N;
      int in = (i+1) % N;
      forces[2*i+0] = K*(midpoints[2*ip+0] + midpoints[2*in+0] - 2*midpoints[2*i+0])*inv_dtheta2;
      forces[2*i+1] = K*(midpoints[2*ip+1] + midpoints[2*in+1] - 2*midpoints[2*i+1])*inv_dtheta2;
    }
  }
private:
  static Array<real_t> circle(int N,
                              real_t radius,
                              real_t center_x,
                              real_t center_y)
  {
    Array<real_t> points(N);
    real_t dtheta = (2. * M_PI) / N;
    for (int i = 0; i < N; i++) {
      points[2 * i + 0] = radius * cos(dtheta * i) + center_x;
      points[2 * i + 1] = radius * sin(dtheta * i) + center_y;
    }
    return points;
  }
};
} // namespace beam