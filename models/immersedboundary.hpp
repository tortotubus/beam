#pragma once

#include "../config/config.hpp"
#include "../general/array.hpp"
#include "../general/error.hpp"

namespace beam {
/**
 * @brief This class provides a Lagrangian mesh for use with the immersed
 * boundary method
 */
class IBStructureMesh
{
protected:
  Array<real_t> points;
  Array<real_t> forces;
  // Array<real_t> velocity;

  int NumberOfPoints, Dimension;

public:
  IBStructureMesh(int Dimension, int NumberOfPoints)
    : Dimension(Dimension)
    , NumberOfPoints(NumberOfPoints)
    , points(Dimension * NumberOfPoints)
    , forces(Dimension * NumberOfPoints)
  //, velocity(Dimension * NumberOfPoints)
  {};

  int GetNumberOfPoints() { return NumberOfPoints; }
  int GetDimension() { return Dimension; }
  Array<real_t>& GetPoints() { return points; }
  Array<real_t>& GetForces() { return forces; }
};

/**
 * @brief This class provides standard methods to help interface with fluid
 * solvers using the immersed boundary method of Peskin.
 */
class IBStructureModel
{
protected:
  IBStructureMesh mesh, mesh_next;

  void CopyCurrentToNext() { mesh_next.GetPoints().Copy(mesh.GetPoints()); }

  virtual void ComputeMidpointForces() = 0;

  virtual void ComputeMidpointPoints(Array<real_t>& velocity, real_t dt)
  {
    int n = mesh_next.GetPoints().Size();
    for (int i = 0; i < n; i++) {
      mesh_next.GetPoints()[i] += 0.5 * dt * velocity[i];
    }
  };

  virtual void ComputeNextPoints(Array<real_t>& velocity, real_t dt)
  {
    Array<real_t> &points = mesh.GetPoints();
    int n = points.Size();
    for (int i = 0; i < n; i++) {
      points[i] += dt * velocity[i];
    }
  };

public:
  IBStructureModel(int Dimension, int NumberOfPoints)
    : mesh(Dimension, NumberOfPoints)
    , mesh_next(Dimension, NumberOfPoints) {};

  int GetNumberOfPoints() { return mesh.GetNumberOfPoints(); }
  int GetDimension() { return mesh.GetDimension(); }

  virtual IBStructureMesh& GetCurrent() {
    return mesh;
  }
  
  virtual IBStructureMesh& GetMidpoint(Array<real_t>& velocity, real_t dt)
  {
    BEAM_ASSERT(velocity.Size() == mesh.GetPoints().Size(),
                "Velocity must match dimension and number of points");

    CopyCurrentToNext();
    ComputeMidpointPoints(velocity, dt);
    ComputeMidpointForces();

    return mesh_next;
  }
  virtual IBStructureMesh& GetNext(Array<real_t>& velocity, real_t dt)
  {
    BEAM_ASSERT(velocity.Size() == mesh.GetPoints().Size(),
                "Velocity must match dimension and number of points");

    ComputeNextPoints(velocity, dt);

    return mesh;
  }
};

} // namespace beam