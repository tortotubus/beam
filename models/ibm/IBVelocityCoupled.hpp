#pragma once

#include "elff/config/config.hpp"
#include "elff/general/error.hpp"

#include "elff/models/ibm/IBMesh.hpp"
#include "elff/models/ibm/IBModel.hpp"

#include <algorithm>
#include <cstdint>
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
class IBVelocityCoupled :
  public IBModel
{
protected:
  IBMesh mesh, mesh_next;

  void CopyCurrentToNext() { mesh_next = mesh; }

  void SetNodalMeasures(const std::vector<real_t>& nodal_measures)
  {
    ELFF_ASSERT(nodal_measures.size() == mesh.GetNumberOfPoints(),
                "Nodal measure length must match number of points.\n");
    std::copy(nodal_measures.begin(), nodal_measures.end(), mesh.GetMeasures().begin());
    std::copy(nodal_measures.begin(),
              nodal_measures.end(),
              mesh_next.GetMeasures().begin());
  }

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

  void pack_state(IBModelState& s) const override
  {
    static constexpr int64_t state_version = 1;

    s.ints.clear();
    s.reals.clear();
    s.bytes.clear();

    const size_t n = mesh.GetNumberOfPoints();
    s.ints.push_back(state_version);
    s.ints.push_back(static_cast<int64_t>(n));

    s.reals.reserve(20 * n);
    pack_mesh_state(s, mesh);
    pack_mesh_state(s, mesh_next);
  }

  void unpack_state(const IBModelState& s) override
  {
    static constexpr int64_t state_version = 1;

    ELFF_VERIFY(s.ints.size() == 2,
                "IBVelocityCoupled::unpack_state(): invalid integer metadata.\n");
    ELFF_VERIFY(s.ints[0] == state_version,
                "IBVelocityCoupled::unpack_state(): unsupported state version.\n");

    const size_t n = mesh.GetNumberOfPoints();
    ELFF_VERIFY(static_cast<size_t>(s.ints[1]) == n,
                "IBVelocityCoupled::unpack_state(): point count mismatch.\n");
    ELFF_VERIFY(s.reals.size() == 20 * n,
                "IBVelocityCoupled::unpack_state(): invalid real buffer size.\n");

    size_t k = 0;
    unpack_mesh_state(s, k, mesh);
    unpack_mesh_state(s, k, mesh_next);
    ELFF_VERIFY(k == s.reals.size(),
                "IBVelocityCoupled::unpack_state(): trailing real data.\n");
  }

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

private:
  static void pack_vertex(IBModelState& s, const IBMesh::IBVertex& value)
  {
    s.reals.push_back(value.x);
    s.reals.push_back(value.y);
    s.reals.push_back(value.z);
  }

  static void unpack_vertex(const IBModelState& s,
                            size_t& k,
                            IBMesh::IBVertex& value)
  {
    value.x = s.reals[k++];
    value.y = s.reals[k++];
    value.z = s.reals[k++];
  }

  static void pack_mesh_state(IBModelState& s, const IBMesh& ib_mesh)
  {
    const size_t n = ib_mesh.GetNumberOfPoints();
    const auto& points = ib_mesh.GetPoints();
    const auto& velocity = ib_mesh.GetVelocity();
    const auto& forces = ib_mesh.GetForces();
    const auto& measures = ib_mesh.GetMeasures();

    for (size_t i = 0; i < n; ++i)
      pack_vertex(s, points[i]);
    for (size_t i = 0; i < n; ++i)
      pack_vertex(s, velocity[i]);
    for (size_t i = 0; i < n; ++i)
      pack_vertex(s, forces[i]);
    for (size_t i = 0; i < n; ++i)
      s.reals.push_back(measures[i]);
  }

  static void unpack_mesh_state(const IBModelState& s,
                                size_t& k,
                                IBMesh& ib_mesh)
  {
    const size_t n = ib_mesh.GetNumberOfPoints();
    auto& points = ib_mesh.GetPoints();
    auto& velocity = ib_mesh.GetVelocity();
    auto& forces = ib_mesh.GetForces();
    auto& measures = ib_mesh.GetMeasures();

    for (size_t i = 0; i < n; ++i)
      unpack_vertex(s, k, points[i]);
    for (size_t i = 0; i < n; ++i)
      unpack_vertex(s, k, velocity[i]);
    for (size_t i = 0; i < n; ++i)
      unpack_vertex(s, k, forces[i]);
    for (size_t i = 0; i < n; ++i)
      measures[i] = s.reals[k++];
  }
};

} // namespace Models
} // namespace ELFF
