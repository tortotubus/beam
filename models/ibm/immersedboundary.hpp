#pragma once

#include "../../config/config.hpp"
#include "../../general/error.hpp"

#include <iostream>
#include <vector>

namespace ELFF {

/**
 * @brief This class provides a Lagrangian mesh for use with the immersed
 * boundary method
 */
class IBStructureMesh
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
  IBStructureMesh(size_t NumberOfPoints)
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
  IBStructureMesh(const IBStructureMesh& other)
    : NumberOfPoints(other.NumberOfPoints)
    , points(other.points)
    , forces(other.forces) {};

  size_t GetNumberOfPoints() { return NumberOfPoints; }
  std::vector<IBVertex>& GetPoints() { return points; }
  std::vector<IBVertex>& GetVelocity() { return velocity; }
  std::vector<IBVertex>& GetForces() { return forces; }
};

/**
 * @brief This class provides standard methods to help interface with fluid
 * solvers using the immersed boundary method of Peskin.
 *
 * This class is designed to be called by a fluid solver and represents some
 * type of structural model immersed in a fluid. In particular, as in the
 * original method of Peskin, we assume a constitutive structural model that
 * takes nodal Lagrangian positions \f(\vec{X}\f) as inputs and returns internal
 * forces \f(\vec{f}\f) at those nodes.
 *
 * At each time step \f(t^n \to t^{n+1}\f), we first evaluate the forces present
 * at time \f(t^{n}\f), e.g.
 *  \f[\vec F^{n}_{i} = \vec F_{\rm internal}(\{\vec X^n_{i}\})\f]
 * for each marker \f(i\f). We then spread the internal constitutive force to
 * the Eulerian grid using a body force density using a discrete delta
 * \f(\delta_h\f) as
 *  \f[\vec f^n(\vec x_k) = \sum_{i} \vec F^n_i \delta_h(\vec x_k - \vec X_i^n)
 * \Delta s_i\f] which is included as an Eulerian forcing term in the momentum
 * equation. Navier-stokes (N-S) is then advanced to the next time step
 * \f(t^{n+1}\f). In the N-S solver, the velocity \f(\vec U^{n+1}_i\f) is then
 * interpolated at the nodes \f(\vec X^n_i\f) using the same \f(\delta_h\f)
 * kernel, e.g.
 *  \f[\vec{U}^{n+1}_i = \sum_{k} \vec u^{n+1}(\vec x_k) \delta_h(\vec x_k -
 * \vec X^n_i)h^d,\f] and the the Lagrangian nodes are then advected according
 * to these velocities. By default, an explicit first-order update is provided
 * as
 *  \f[\vec X^{n+1}_i = \vec X^n_i + \Delta t \vec{U}^{n+1}_i.\f]
 */
class IBVelocityCoupledStructureModel
{
protected:
  IBStructureMesh mesh, mesh_next;

  void CopyCurrentToNext() { mesh_next = mesh; }

  virtual void ComputeMidpointForces() = 0;

  virtual void ComputeMidpointPoints(
    std::vector<IBStructureMesh::IBVertex>& velocity,
    real_t dt)
  {
    std::vector<IBStructureMesh::IBVertex>& midpoints = mesh_next.GetPoints();
    ELFF_ASSERT(velocity.size() == midpoints.size(),
                "Velocity vector must have same size as mesh.\n");

    for (size_t pi = 0; pi < midpoints.size(); pi++) {
      midpoints[pi].x += 0.5 * dt * velocity[pi].x;
      midpoints[pi].y += 0.5 * dt * velocity[pi].y;
      midpoints[pi].z += 0.5 * dt * velocity[pi].z;
    }
  };

  virtual void ComputeNextPoints(
    std::vector<IBStructureMesh::IBVertex>& velocity,
    real_t dt)
  {
    std::vector<IBStructureMesh::IBVertex>& midpoints = mesh.GetPoints();
    ELFF_ASSERT(velocity.size() == midpoints.size(),
                "Velocity vector must have same size as mesh.\n");

    for (size_t pi = 0; pi < midpoints.size(); pi++) {
      midpoints[pi].x += dt * velocity[pi].x;
      midpoints[pi].y += dt * velocity[pi].y;
      midpoints[pi].z += dt * velocity[pi].z;
    }
  };

public:
  IBVelocityCoupledStructureModel(size_t NumberOfPoints)
    : mesh(NumberOfPoints)
    , mesh_next(NumberOfPoints) {};

  size_t GetNumberOfPoints() { return mesh.GetNumberOfPoints(); }

  virtual IBStructureMesh& GetCurrent() { return mesh; }

  virtual IBStructureMesh& GetMidpoint(
    std::vector<IBStructureMesh::IBVertex>& velocity,
    real_t dt)
  {
    ELFF_ASSERT(velocity.size() == mesh.GetPoints().size(),
                "Velocity must match dimension and number of points");

    CopyCurrentToNext();
    ComputeMidpointPoints(velocity, dt);
    ComputeMidpointForces();

    return mesh_next;
  }

  virtual IBStructureMesh& GetMidpoint(IBStructureMesh::IBVertex* velocity,
                                       size_t n,
                                       real_t dt)
  {
    std::vector<IBStructureMesh::IBVertex> velocity_wrapper(velocity,
                                                            velocity + n);
    return GetMidpoint(velocity_wrapper, dt);
  }

  virtual IBStructureMesh& GetNext(
    std::vector<IBStructureMesh::IBVertex>& velocity,
    real_t dt)
  {
    ELFF_ASSERT(velocity.size() == mesh.GetPoints().size(),
                "Velocity must match dimension and number of points");

    ComputeNextPoints(velocity, dt);

    return mesh;
  }

  virtual IBStructureMesh& GetNext(IBStructureMesh::IBVertex* velocity,
                                   size_t n,
                                   real_t dt)
  {
    std::vector<IBStructureMesh::IBVertex> velocity_wrapper(velocity,
                                                            velocity + n);
    return GetNext(velocity_wrapper, dt);
  }
};

/**
 * @brief Immersed boundary coupling abstract base class for constitutive models
 * that take forces as inputs and give position and velocity as outputs
 */
class IBForceCoupledStructureModel
{
protected:
  IBStructureMesh mesh, mesh_next;

  void CopyCurrentToNext() { mesh_next = mesh; }

  virtual void ComputeNextPoints(std::vector<IBStructureMesh::IBVertex> force,
                                 real_t dt) = 0;

public:
  IBForceCoupledStructureModel(size_t NumberOfPoints)
    : mesh(NumberOfPoints)
    , mesh_next(NumberOfPoints) {};

  size_t GetNumberOfPoints()
  {
    // std::cout << "IBForceCoupledStructureModel::GetNumberOfPoint()\n";
    return mesh.GetNumberOfPoints();
  }

  virtual IBStructureMesh& GetCurrent() { return mesh; }

  virtual IBStructureMesh& GetNext(IBStructureMesh::IBVertex* force,
                                   size_t n,
                                   real_t dt)
  {
    std::vector<IBStructureMesh::IBVertex> force_wrapper(force, force + n);
    return GetNext(force_wrapper, dt);
  }

  virtual IBStructureMesh& GetNext(std::vector<IBStructureMesh::IBVertex> force,
                                   real_t dt)
  {
    ELFF_ASSERT(force.size() == mesh.GetPoints().size(),
                "Force length must match number of points.\n");

    ComputeNextPoints(force, dt);

    return mesh_next;
  }
};

} // namespace ELFF