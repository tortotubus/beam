#include "elff/models/beam/IBEulerBeamPenalty.hpp"

namespace ELFF {
namespace Models {

void
IBEulerBeamPenalty::EBMeshToIBMeshNext()
{
  EulerBeamInextensiblePenalty::update_mesh();

  auto& ib_points = IBForceCoupled::mesh_next.GetPoints();
  auto& ib_velocity = IBForceCoupled::mesh_next.GetVelocity();

  auto& eb_points = EulerBeam::mesh.get_centerline();
  auto& eb_velocity = EulerBeam::mesh.get_centerline_velocity();

  for (size_t ni = 0; ni < nodes; ++ni) {
    ib_points[ni].x = eb_points[ni][0];
    ib_points[ni].y = eb_points[ni][1];
    ib_points[ni].z = eb_points[ni][2];

    ib_velocity[ni].x = eb_velocity[ni][0];
    ib_velocity[ni].y = eb_velocity[ni][1];
    ib_velocity[ni].z = eb_velocity[ni][2];
  }
}

void
IBEulerBeamPenalty::EBMeshToIBMeshCurrent()
{
  EulerBeamInextensiblePenalty::update_mesh();

  auto& ib_points = IBForceCoupled::mesh.GetPoints();
  auto& ib_velocity = IBForceCoupled::mesh.GetVelocity();

  auto& eb_points = EulerBeam::mesh.get_centerline();
  auto& eb_velocity = EulerBeam::mesh.get_centerline_velocity();

  for (size_t ni = 0; ni < nodes; ++ni) {
    ib_points[ni].x = eb_points[ni][0];
    ib_points[ni].y = eb_points[ni][1];
    ib_points[ni].z = eb_points[ni][2];

    ib_velocity[ni].x = eb_velocity[ni][0];
    ib_velocity[ni].y = eb_velocity[ni][1];
    ib_velocity[ni].z = eb_velocity[ni][2];
  }
}

IBEulerBeamPenalty::IBEulerBeamPenalty(real_t length,
                                       real_t EI,
                                       real_t mu,
                                       size_t nodes,
                                       EulerBeamBCs bcs,
                                       real_t r_penalty)
  : EulerBeamInextensiblePenalty(length, EI, mu, nodes, bcs, r_penalty)
  , IBForceCoupled(nodes)
{
  const size_t beam_nodes = EulerBeam::mesh.get_nodes();
  const real_t ds = EulerBeam::mesh.get_ds();
  std::vector<real_t> measures(beam_nodes, ds);
  if (beam_nodes > 1) {
    measures.front() *= 0.5;
    measures.back() *= 0.5;
  }
  SetNodalMeasures(measures);
  EBMeshToIBMeshCurrent();
  EBMeshToIBMeshNext();
}

void
IBEulerBeamPenalty::apply_initial_condition(EulerBeamMesh& mesh)
{
  EulerBeamInextensiblePenalty::apply_initial_condition(mesh);
  EBMeshToIBMeshCurrent();
  EBMeshToIBMeshNext();
}

void
IBEulerBeamPenalty::ComputeNextPoints(std::vector<IBMesh::IBVertex> force,
                                      real_t                       dt)
{
  ELFF_ASSERT(force.size() == nodes, "Force array must match node count.\n");

  std::vector<std::array<real_t, 3>> load(nodes);
  for (size_t ni = 0; ni < nodes; ++ni) {
    load[ni][0] = force[ni].x;
    load[ni][1] = force[ni].y;
    load[ni][2] = force[ni].z;
  }

  EulerBeamInextensiblePenalty::solve(dt, load);
  EBMeshToIBMeshNext();
}

void
IBEulerBeamPenalty::pack_state(IBModelState& s) const
{
  static constexpr int64_t state_version = 1;

  s.ints.clear();
  s.reals.clear();
  s.bytes.clear();

  s.ints.push_back(state_version);
  s.ints.push_back(static_cast<int64_t>(nodes));
  s.ints.push_back(static_cast<int64_t>(ndof));
  s.ints.push_back(static_cast<int64_t>(time_iter));

  s.reals.reserve(1 + 3 * ndof);
  s.reals.push_back(t);

  auto pack_vec = [&](const VectorXd& v) {
    for (Eigen::Index i = 0; i < v.size(); ++i)
      s.reals.push_back(v(i));
  };

  pack_vec(u);
  pack_vec(v_prev);
  pack_vec(a_prev);
}

void
IBEulerBeamPenalty::unpack_state(const IBModelState& s)
{
  static constexpr int64_t state_version = 1;

  ELFF_ASSERT(s.ints.size() == 4,
              "IBEulerBeamPenalty::unpack_state(): invalid integer metadata.\n");
  ELFF_ASSERT(s.ints[0] == state_version,
              "IBEulerBeamPenalty::unpack_state(): unsupported state version.\n");
  ELFF_ASSERT(static_cast<size_t>(s.ints[1]) == nodes,
              "IBEulerBeamPenalty::unpack_state(): node count mismatch.\n");
  ELFF_ASSERT(static_cast<size_t>(s.ints[2]) == ndof,
              "IBEulerBeamPenalty::unpack_state(): ndof mismatch.\n");

  const size_t expected_reals = 1 + 3 * ndof;
  ELFF_ASSERT(s.reals.size() == expected_reals,
              "IBEulerBeamPenalty::unpack_state(): invalid real buffer size.\n");

  size_t k = 0;
  t = s.reals[k++];

  auto unpack_vec = [&](VectorXd& v) {
    for (Eigen::Index i = 0; i < v.size(); ++i)
      v(i) = s.reals[k++];
  };

  unpack_vec(u);
  unpack_vec(v_prev);
  unpack_vec(a_prev);

  time_iter = static_cast<size_t>(s.ints[3]);
  u_prev = u;

  EulerBeamInextensiblePenalty::update_mesh();
  EBMeshToIBMeshCurrent();
  EBMeshToIBMeshNext();
}

} // namespace Models
} // namespace ELFF
