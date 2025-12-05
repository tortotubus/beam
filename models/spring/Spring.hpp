#pragma once

#include "../ibm/immersedboundary.hpp"

#include <cmath>

namespace beam {
class IBSpringCircle : public IBVelocityCoupledStructureModel
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
    : IBVelocityCoupledStructureModel(NumberOfPoints)
    , K(K)
    , radius_x(radius)
    , radius_y(radius)
    , dtheta((2 * M_1_PI) / NumberOfPoints)
  {
    std::vector<IBStructureMesh::IBVertex> circle_points =
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
    : IBVelocityCoupledStructureModel(NumberOfPoints)
    , K(K)
    , radius_x(radius_x)
    , radius_y(radius_y)
    , dtheta((2 * M_1_PI) / NumberOfPoints)
  {
    std::vector<IBStructureMesh::IBVertex> ellipse_point =
      ellipse(NumberOfPoints, radius_x, radius_y, center_x, center_y);

    for (int i = 0; i < NumberOfPoints; ++i) {
      mesh.GetPoints()[i] = ellipse_point[i];
    }
  };

protected:
  void ComputeMidpointForces() override
  {
    // Array<fem::Vertex>& forces = mesh_next.GetForces();
    // Array<fem::Vertex>& midpoints = mesh_next.GetPoints();

    std::vector<IBStructureMesh::IBVertex>& forces = mesh_next.GetForces();
    std::vector<IBStructureMesh::IBVertex>& midpoints = mesh_next.GetPoints();

    int N = GetNumberOfPoints();
    real_t inv_dtheta2 = 1. / (dtheta * dtheta);

    for (int i = 0; i < N; i++) {
      int ip = (i - 1 + N) % N;
      int in = (i + 1) % N;
      forces[i].x = K *
                    (midpoints[ip].x + midpoints[in].x - 2 * midpoints[i].x) *
                    inv_dtheta2;
      forces[i].y = K *
                    (midpoints[ip].y + midpoints[in].y - 2 * midpoints[i].y) *
                    inv_dtheta2;
    }
  }

private:
  static std::vector<IBStructureMesh::IBVertex> circle(int N,
                                                       real_t radius,
                                                       real_t center_x,
                                                       real_t center_y)
  {
    std::vector<IBStructureMesh::IBVertex> points(N);
    // Array<fem::Vertex> points(N);
    real_t dtheta = (2. * M_PI) / N;
    for (int pi = 0; pi < N; pi++) {
      points[pi].x = radius * cos(dtheta * pi) + center_x;
      points[pi].y = radius * sin(dtheta * pi) + center_y;
    }
    return points;
  }

  static std::vector<IBStructureMesh::IBVertex> ellipse(int N,
                                                        real_t radius_x,
                                                        real_t radius_y,
                                                        real_t center_x,
                                                        real_t center_y)
  {
    // Array<fem::Vertex> points(N);
    std::vector<IBStructureMesh::IBVertex> points(N);
    real_t dtheta = (2. * M_PI) / N;
    for (int pi = 0; pi < N; pi++) {
      points[pi].x = radius_x * cos(dtheta * pi) + center_x;
      points[pi].y = radius_y * sin(dtheta * pi) + center_y;
    }
    return points;
  }
};
} // namespace beam