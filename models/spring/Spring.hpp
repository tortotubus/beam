#pragma once

#include "../ibm/immersedboundary.hpp"

#include <cmath>

namespace beam {

/**
 * @class IBSpringCircle
 * @brief A quasi-static loop of springs arranged in a circular or elliptical
 * configuration.
 *
 * This class models a structure composed of N points arranged in a circle (or
 * ellipse) connected by linear springs. The springs apply elastic restoring
 * forces based on the discrete Laplacian of the point positions, creating a
 * quasi-static equilibrium when coupled with fluid velocity.
 *
 * Each spring resists changes in distance between adjacent points. The force
 * on point i depends on its position relative to its neighbors (i-1) and
 * (i+1), scaled by the spring stiffness K and the angular spacing dθ.
 *
 * The model supports both circular (uniform radius) and elliptical
 * (different radii in x and y) initial configurations.
 *
 * @note This class inherits from IBVelocityCoupledStructureModel and
 *       implements a velocity-coupled immersed boundary method compatible
 *       with fluid solvers.
 */
class IBSpringCircle : public IBVelocityCoupledStructureModel
{
protected:
  real_t K;              ///< Spring stiffness coefficient
  real_t radius_x;       ///< Radius (or semi-axis) in x-direction
  real_t radius_y;       ///< Radius (or semi-axis) in y-direction
  real_t dtheta;         ///< Angular spacing between adjacent points (2π / N)

public:
  /**
   * @brief Constructor for a circular spring loop.
   *
   * Creates a quasi-static spring loop with N points equally spaced on a
   * circle of uniform radius. All points are connected to their neighbors
   * via linear springs with stiffness K.
   *
   * @param NumberOfPoints Number of spring masses in the loop
   * @param K              Spring stiffness coefficient
   * @param radius         Radius of the circular configuration
   * @param center_x       x-coordinate of the circle center
   * @param center_y       y-coordinate of the circle center
   */
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

  /**
   * @brief Constructor for an elliptical spring loop.
   *
   * Creates a quasi-static spring loop with N points equally spaced on an
   * ellipse with potentially different semi-axes in the x and y directions.
   * All points are connected to their neighbors via linear springs with
   * stiffness K.
   *
   * @param NumberOfPoints Number of spring masses in the loop
   * @param K              Spring stiffness coefficient
   * @param radius_x       Semi-axis length in the x-direction
   * @param radius_y       Semi-axis length in the y-direction
   * @param center_x       x-coordinate of the ellipse center
   * @param center_y       y-coordinate of the ellipse center
   */
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
  /**
   * @brief Compute elastic spring forces at each midpoint.
   *
   * Evaluates the discrete Laplacian-based spring forces on all points in
   * the structure. The force on point i is proportional to the curvature
   * (second difference) of the shape:
   *
   *     F_i = K * (x_{i-1} + x_{i+1} - 2*x_i) / (dθ)²
   *
   * This assumes a quasi-static configuration where springs resist local
   * deformations from the rest shape. Periodic boundary conditions connect
   * the last point to the first.
   *
   * @see IBVelocityCoupledStructureModel::ComputeMidpointForces()
   */
  void ComputeMidpointForces() override
  {
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
  /**
   * @brief Generate N points on a circle.
   *
   * Creates an evenly spaced set of points around a circle of given radius,
   * centered at (center_x, center_y).
   *
   * @param N       Number of points to generate
   * @param radius  Radius of the circle
   * @param center_x x-coordinate of the center
   * @param center_y y-coordinate of the center
   *
   * @return A vector of N IBVertex points on the circle
   */
  static std::vector<IBStructureMesh::IBVertex> circle(int N,
                                                       real_t radius,
                                                       real_t center_x,
                                                       real_t center_y)
  {
    std::vector<IBStructureMesh::IBVertex> points(N);
    real_t dtheta = (2. * M_PI) / N;
    for (int pi = 0; pi < N; pi++) {
      points[pi].x = radius * cos(dtheta * pi) + center_x;
      points[pi].y = radius * sin(dtheta * pi) + center_y;
    }
    return points;
  }

  /**
   * @brief Generate N points on an ellipse.
   *
   * Creates an evenly spaced set of points around an ellipse with semi-axes
   * radius_x and radius_y, centered at (center_x, center_y).
   *
   * @param N        Number of points to generate
   * @param radius_x Semi-axis length in the x-direction
   * @param radius_y Semi-axis length in the y-direction
   * @param center_x x-coordinate of the center
   * @param center_y y-coordinate of the center
   *
   * @return A vector of N IBVertex points on the ellipse
   */
  static std::vector<IBStructureMesh::IBVertex> ellipse(int N,
                                                        real_t radius_x,
                                                        real_t radius_y,
                                                        real_t center_x,
                                                        real_t center_y)
  {
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