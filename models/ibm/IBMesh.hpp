#pragma once

#include <cstddef>
#include <vector>

#include "elff/config/config.hpp"

namespace ELFF {
namespace Models {

/**
 * @brief This class provides a Lagrangian mesh for use with the immersed
 * boundary method
 */
class IBMesh
{
public:
  /**
   * @brief One immersed-boundary mesh vertex stored as a 3D Cartesian point.
   *
   * The same lightweight coordinate container is reused for IB positions,
   * velocities, and forces throughout the mesh interface.
   */
  struct IBVertex
  {
    real_t x, y, z; ///< Cartesian vector components.
  };

protected:
  size_t NumberOfPoints;

  std::vector<IBVertex> points;
  std::vector<IBVertex> forces;
  std::vector<IBVertex> velocity;
  std::vector<real_t> measures;

public:
  IBMesh(size_t NumberOfPoints)
    : NumberOfPoints(NumberOfPoints)
    , points(NumberOfPoints, { 0, 0, 0 })
    , forces(NumberOfPoints, { 0, 0, 0 })
    , velocity(NumberOfPoints, { 0, 0, 0 })
    , measures(NumberOfPoints, 1.0) {};

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
    , forces(other.forces)
    , velocity(other.velocity)
    , measures(other.measures) {};

  size_t GetNumberOfPoints() { return NumberOfPoints; }
  size_t GetNumberOfPoints() const { return NumberOfPoints; }
  std::vector<IBVertex>& GetPoints() { return points; }
  const std::vector<IBVertex>& GetPoints() const { return points; }
  std::vector<IBVertex>& GetVelocity() { return velocity; }
  const std::vector<IBVertex>& GetVelocity() const { return velocity; }
  std::vector<IBVertex>& GetForces() { return forces; }
  const std::vector<IBVertex>& GetForces() const { return forces; }
  std::vector<real_t>& GetMeasures() { return measures; }
  const std::vector<real_t>& GetMeasures() const { return measures; }
};


} // namespace Models
} // namespace ELFF
