#pragma once

#include <vector>
#include "config/config.hpp"

namespace ELFF {
namespace Models {

/**
 * @brief This class provides a Lagrangian mesh for use with the immersed
 * boundary method
 */
class IBMesh
{
public:
  struct IBVertex
  {
    real_t x, y, z;
  };

protected:
  size_t NumberOfPoints;

  std::vector<IBVertex> points;
  std::vector<IBVertex> forces;
  std::vector<IBVertex> velocity;

public:
  IBMesh(size_t NumberOfPoints)
    : NumberOfPoints(NumberOfPoints)
    , points(NumberOfPoints, { 0, 0, 0 })
    , forces(NumberOfPoints, { 0, 0, 0 })
    , velocity(NumberOfPoints, { 0, 0, 0 }) {};

  /**
   * @brief Copy constructor
   *
   * Performs a memberwise copy of the mesh data. This creates a separate
   * container for points and forces so the copy can be modified without
   * affecting the original object.
   */
  IBMesh(const IBMesh& other)
    : NumberOfPoints(other.NumberOfPoints)
    , points(other.points)
    , forces(other.forces) {};

  size_t GetNumberOfPoints() { return NumberOfPoints; }
  std::vector<IBVertex>& GetPoints() { return points; }
  std::vector<IBVertex>& GetVelocity() { return velocity; }
  std::vector<IBVertex>& GetForces() { return forces; }
};


} // namespace Models
} // namespace ELFF
