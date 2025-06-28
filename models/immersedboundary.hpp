#pragma once

#include "../config/config.hpp"
#include "../general/array.hpp"
#include "../mesh/vertex.hpp"
#include "../general/error.hpp"

namespace beam {
/**
 * @brief This class provides a Lagrangian mesh for use with the immersed
 * boundary method
 */
class IBStructureMesh
{
protected:
  Array<fem::Vertex> points;
  Array<fem::Vertex> forces;

  int NumberOfPoints, Dimension;

public:
  IBStructureMesh(int Dimension, int NumberOfPoints)
    : Dimension(Dimension)
    , NumberOfPoints(NumberOfPoints)
    , points(NumberOfPoints)
    , forces(NumberOfPoints)
  //, velocity(Dimension * NumberOfPoints)
  {};

  int GetNumberOfPoints() { return NumberOfPoints; }
  int GetDimension() { return Dimension; }
  Array<fem::Vertex>& GetPoints() { return points; }
  Array<fem::Vertex>& GetForces() { return forces; }
};

/**
 * @brief This class provides standard methods to help interface with fluid
 * solvers using the immersed boundary method of Peskin.
 */
class IBStructureModel
{
protected:
  IBStructureMesh mesh, mesh_next;

  void CopyCurrentToNext() { 
    //mesh_next.GetPoints().Copy(mesh.GetPoints()); 
    size_t nn = mesh.GetNumberOfPoints();
    for(size_t i = 0; i < nn; ++i) {
      mesh_next.GetPoints()[i](0) = mesh.GetPoints()[i](0);
      mesh_next.GetPoints()[i](1) = mesh.GetPoints()[i](1);
      mesh_next.GetPoints()[i](2) = mesh.GetPoints()[i](2);
    }
  }

  virtual void ComputeMidpointForces() = 0;

  virtual void ComputeMidpointPoints(Array<fem::Vertex>& velocity, real_t dt)
  {
    Array<fem::Vertex> &midpoints = mesh_next.GetPoints();
    int n = midpoints.Size();
    int d = mesh_next.GetDimension();

    for (int pi = 0; pi < n; pi++) {
      midpoints[pi](0) += 0.5 * dt * velocity[pi](0);
      midpoints[pi](1) += 0.5 * dt * velocity[pi](1);
      midpoints[pi](2) += 0.5 * dt * velocity[pi](2);
    }
  };

  virtual void ComputeNextPoints(Array<fem::Vertex>& velocity, real_t dt)
  {
    Array<fem::Vertex> &points = mesh.GetPoints();
    int n = points.Size();
    int d = mesh.GetDimension();
    for (int pi = 0; pi < n; pi++) {
      points[pi](0) += dt * velocity[pi](0);
      points[pi](1) += dt * velocity[pi](1);
      points[pi](2) += dt * velocity[pi](2);
    }
  };

public:
  IBStructureModel(int Dimension, int NumberOfPoints)
    : mesh(Dimension, NumberOfPoints)
    , mesh_next(Dimension, NumberOfPoints) {};

  int GetNumberOfPoints() { return mesh.GetNumberOfPoints(); }
  int GetDimension() { return mesh.GetDimension(); }

  virtual IBStructureMesh& GetCurrent() { return mesh; }
  
  virtual IBStructureMesh& GetMidpoint(Array<fem::Vertex>& velocity, real_t dt)
  {
    BEAM_ASSERT(velocity.Size() == mesh.GetPoints().Size(),
                "Velocity must match dimension and number of points");

    CopyCurrentToNext();
    ComputeMidpointPoints(velocity, dt);
    ComputeMidpointForces();

    return mesh_next;
  }

  virtual IBStructureMesh& GetMidpoint(fem::Vertex *velocity, size_t n, real_t dt) {
    Array<fem::Vertex> velocity_wrapper(velocity, n, false);
    return GetMidpoint(velocity_wrapper, dt);
  }

  virtual IBStructureMesh& GetNext(Array<fem::Vertex>& velocity, real_t dt)
  {
    BEAM_ASSERT(velocity.Size() == mesh.GetPoints().Size(),
                "Velocity must match dimension and number of points");

    ComputeNextPoints(velocity, dt);

    return mesh;
  }

  virtual IBStructureMesh& GetNext(fem::Vertex *velocity, size_t n, real_t dt) {
    Array<fem::Vertex> velocity_wrapper(velocity, n, false);
    return GetNext(velocity_wrapper, dt);
  }

};

} // namespace beam