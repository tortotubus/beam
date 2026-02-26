#pragma once

#include "config/config.hpp"
#include "general/error.hpp"

#include "models/ibm/IBMesh.hpp"

#include <vector>

namespace ELFF {
namespace Models {

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
class IBVelocityCoupled
{
protected:
  IBMesh mesh, mesh_next;

  void CopyCurrentToNext() { mesh_next = mesh; }

  virtual void ComputeMidpointForces() = 0;

  virtual void ComputeMidpointPoints(
    std::vector<IBMesh::IBVertex>& velocity,
    real_t dt)
  {
    std::vector<IBMesh::IBVertex>& midpoints = mesh_next.GetPoints();
    ELFF_ASSERT(velocity.size() == midpoints.size(),
                "Velocity vector must have same size as mesh.\n");

    for (size_t pi = 0; pi < midpoints.size(); pi++) {
      midpoints[pi].x += 0.5 * dt * velocity[pi].x;
      midpoints[pi].y += 0.5 * dt * velocity[pi].y;
      midpoints[pi].z += 0.5 * dt * velocity[pi].z;
    }
  };

  virtual void ComputeNextPoints(
    std::vector<IBMesh::IBVertex>& velocity,
    real_t dt)
  {
    std::vector<IBMesh::IBVertex>& midpoints = mesh.GetPoints();
    ELFF_ASSERT(velocity.size() == midpoints.size(),
                "Velocity vector must have same size as mesh.\n");

    for (size_t pi = 0; pi < midpoints.size(); pi++) {
      midpoints[pi].x += dt * velocity[pi].x;
      midpoints[pi].y += dt * velocity[pi].y;
      midpoints[pi].z += dt * velocity[pi].z;
    }
  };

public:
  IBVelocityCoupled(size_t NumberOfPoints)
    : mesh(NumberOfPoints)
    , mesh_next(NumberOfPoints) {};

  size_t GetNumberOfPoints() { return mesh.GetNumberOfPoints(); }

  virtual IBMesh& GetCurrent() { return mesh; }

  virtual IBMesh& GetMidpoint(
    std::vector<IBMesh::IBVertex>& velocity,
    real_t dt)
  {
    ELFF_ASSERT(velocity.size() == mesh.GetPoints().size(),
                "Velocity must match dimension and number of points");

    CopyCurrentToNext();
    ComputeMidpointPoints(velocity, dt);
    ComputeMidpointForces();

    return mesh_next;
  }

  virtual IBMesh& GetMidpoint(IBMesh::IBVertex* velocity,
                                       size_t n,
                                       real_t dt)
  {
    std::vector<IBMesh::IBVertex> velocity_wrapper(velocity,
                                                            velocity + n);
    return GetMidpoint(velocity_wrapper, dt);
  }

  virtual IBMesh& GetNext(
    std::vector<IBMesh::IBVertex>& velocity,
    real_t dt)
  {
    ELFF_ASSERT(velocity.size() == mesh.GetPoints().size(),
                "Velocity must match dimension and number of points");

    ComputeNextPoints(velocity, dt);

    return mesh;
  }

  virtual IBMesh& GetNext(IBMesh::IBVertex* velocity,
                                   size_t n,
                                   real_t dt)
  {
    std::vector<IBMesh::IBVertex> velocity_wrapper(velocity,
                                                            velocity + n);
    return GetNext(velocity_wrapper, dt);
  }
};

} // namespace Models
} // namespace ELFF