#include "elff/models/articulated/IBBendingRodImplicit.hpp"
#include <algorithm>
#include <cmath>
#include <climits>
#include "elff/general/error.hpp"

namespace ELFF::Models {
namespace {
size_t nodeCount(int links) {
  if (links < 1 || links > INT_MAX / 6 - 1)
    ELFF_ABORT("invalid rod link count");
  return static_cast<size_t>(links) + 1;
}
}

IBBendingRodImplicit::IBBendingRodImplicit(
    const BendingRodImplicitConfig &config)
  : IBForceCoupled(nodeCount(config.links)), BendingRodImplicit(config)
{
  std::vector<real_t> measures(nodeCount(config.links), get_mesh().get_ds());
  measures.front() *= 0.5;
  measures.back() *= 0.5;
  SetNodalMeasures(measures);
  syncMeshes();
}

void IBBendingRodImplicit::syncMeshes()
{
  const auto &rod_mesh = get_mesh();
  for (size_t i = 0; i < mesh.GetNumberOfPoints(); ++i) {
    const auto &x = rod_mesh.get_centerline(i);
    const auto &v = rod_mesh.get_centerline_velocity(i);
    mesh.GetPoints()[i] = {x[0], x[1], x[2]};
    mesh.GetVelocity()[i] = {v[0], v[1], v[2]};
  }
  CopyCurrentToNext();
}

void IBBendingRodImplicit::setState(const Eigen::VectorXd &state)
{
  BendingRodImplicit::setState(state);
  std::fill(mesh.GetForces().begin(), mesh.GetForces().end(), IBMesh::IBVertex{0, 0, 0});
  syncMeshes();
}

void IBBendingRodImplicit::ComputeNextPoints(
    std::vector<IBMesh::IBVertex> force, real_t dt)
{
  if (force.size() != mesh.GetNumberOfPoints())
    ELFF_ABORT("force count must match rod endpoints");
  std::vector<Vector> density(force.size());
  for (size_t i = 0; i < force.size(); ++i)
    density[i] = {force[i].x, force[i].y, force[i].z};
  solve(dt, density);
}

void IBBendingRodImplicit::solve(real_t dt, const std::vector<Vector>& density)
{
  BendingRodImplicit::solve(dt, density);
  for (size_t i = 0; i < density.size(); ++i)
    mesh.GetForces()[i] = {density[i][0], density[i][1], density[i][2]};
  syncMeshes();
}

void IBBendingRodImplicit::pack_state(IBModelState &s) const
{
  BendingRodImplicit::pack_state(s);
  s.ints.push_back(1); // IB adapter payload version.
  for (const auto &f : mesh.GetForces()) {
    s.reals.push_back(f.x);
    s.reals.push_back(f.y);
    s.reals.push_back(f.z);
  }
}

void IBBendingRodImplicit::unpack_state(const IBModelState &s)
{
  const size_t count = 3 * mesh.GetNumberOfPoints();
  if (s.ints.size() != 6 || s.ints.back() != 1 || s.reals.size() < count)
    ELFF_ABORT("incompatible implicit rod IB checkpoint");
  for (auto value : s.reals)
    if (!std::isfinite(value)) ELFF_ABORT("non-finite implicit rod IB checkpoint");
  IBModelState structural = s;
  structural.ints.pop_back();
  structural.reals.resize(s.reals.size() - count);
  BendingRodImplicit::unpack_state(structural);
  size_t k = structural.reals.size();
  for (auto &f : mesh.GetForces()) {
    f = {s.reals[k], s.reals[k + 1], s.reals[k + 2]};
    k += 3;
  }
  syncMeshes();
}
} // namespace ELFF::Models
