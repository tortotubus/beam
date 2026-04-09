#include "elff/models/beam/IBEulerBeamHuang.hpp"

namespace ELFF {
namespace Models {

namespace {
void
pack_matrix(IBModelState& s, const EulerBeamInextensibleHuang::MatX3& m)
{
  for (Eigen::Index i = 0; i < m.rows(); ++i) {
    for (Eigen::Index j = 0; j < m.cols(); ++j) {
      s.reals.push_back(m(i, j));
    }
  }
}

void
unpack_matrix(const IBModelState&                  s,
              size_t&                              k,
              EulerBeamInextensibleHuang::MatX3&   m)
{
  for (Eigen::Index i = 0; i < m.rows(); ++i) {
    for (Eigen::Index j = 0; j < m.cols(); ++j) {
      m(i, j) = s.reals[k++];
    }
  }
}
} // namespace

void
IBEulerBeamHuang::EBMeshToIBMeshNext()
{
  EulerBeamInextensibleHuang::update_mesh();

  auto& ib_points = IBForceCoupled::mesh_next.GetPoints();
  auto& ib_velocity = IBForceCoupled::mesh_next.GetVelocity();

  auto& eb_points = EulerBeam::mesh.get_centerline();
  auto& eb_velocity = EulerBeam::mesh.get_centerline_velocity();
  const size_t beam_nodes = EulerBeam::mesh.get_nodes();

  for (size_t ni = 0; ni < beam_nodes; ++ni) {
    ib_points[ni].x = eb_points[ni][0];
    ib_points[ni].y = eb_points[ni][1];
    ib_points[ni].z = eb_points[ni][2];

    ib_velocity[ni].x = eb_velocity[ni][0];
    ib_velocity[ni].y = eb_velocity[ni][1];
    ib_velocity[ni].z = eb_velocity[ni][2];
  }
}

void
IBEulerBeamHuang::EBMeshToIBMeshCurrent()
{
  EulerBeamInextensibleHuang::update_mesh();

  auto& ib_points = IBForceCoupled::mesh.GetPoints();
  auto& ib_velocity = IBForceCoupled::mesh.GetVelocity();

  auto& eb_points = EulerBeam::mesh.get_centerline();
  auto& eb_velocity = EulerBeam::mesh.get_centerline_velocity();
  const size_t beam_nodes = EulerBeam::mesh.get_nodes();

  for (size_t ni = 0; ni < beam_nodes; ++ni) {
    ib_points[ni].x = eb_points[ni][0];
    ib_points[ni].y = eb_points[ni][1];
    ib_points[ni].z = eb_points[ni][2];

    ib_velocity[ni].x = eb_velocity[ni][0];
    ib_velocity[ni].y = eb_velocity[ni][1];
    ib_velocity[ni].z = eb_velocity[ni][2];
  }
}

IBEulerBeamHuang::IBEulerBeamHuang(real_t length,
                                   real_t EI,
                                   real_t mu,
                                   size_t nodes,
                                   EulerBeamBCs bcs)
  : EulerBeamInextensibleHuang(length, EI, mu, nodes, bcs)
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
IBEulerBeamHuang::apply_initial_condition(EulerBeamMesh& mesh)
{
  EulerBeamInextensibleHuang::apply_initial_condition(mesh);
  EBMeshToIBMeshCurrent();
  EBMeshToIBMeshNext();
}

void
IBEulerBeamHuang::ComputeNextPoints(std::vector<IBMesh::IBVertex> force,
                                    real_t                       dt)
{
  const size_t beam_nodes = EulerBeam::mesh.get_nodes();
  ELFF_ASSERT(force.size() == beam_nodes,
              "Force array must match node count.\n");

  std::vector<std::array<real_t, 3>> load(beam_nodes);
  for (size_t ni = 0; ni < beam_nodes; ++ni) {
    load[ni][0] = force[ni].x;
    load[ni][1] = force[ni].y;
    load[ni][2] = force[ni].z;
  }

  EulerBeamInextensibleHuang::solve(dt, load);
  EBMeshToIBMeshNext();
}

void
IBEulerBeamHuang::pack_state(IBModelState& s) const
{
  static constexpr int64_t state_version = 1;
  const int64_t beam_nodes = static_cast<int64_t>(EulerBeam::mesh.get_nodes());

  s.ints.clear();
  s.reals.clear();
  s.bytes.clear();

  s.ints.push_back(state_version);
  s.ints.push_back(beam_nodes);
  s.ints.push_back(static_cast<int64_t>(time_iter));
  s.ints.push_back(static_cast<int64_t>(initial_velocity_pending));
  s.ints.push_back(static_cast<int64_t>(implicit_bending));

  s.reals.reserve(2 + 12 * static_cast<size_t>(beam_nodes) +
                  static_cast<size_t>(last_T.size()));
  s.reals.push_back(t);
  s.reals.push_back(last_dt);
  pack_matrix(s, X_init);
  pack_matrix(s, X_n);
  pack_matrix(s, X_nm1);
  pack_matrix(s, V_init);
  for (Eigen::Index i = 0; i < last_T.size(); ++i) {
    s.reals.push_back(last_T(i));
  }
}

void
IBEulerBeamHuang::unpack_state(const IBModelState& s)
{
  static constexpr int64_t state_version = 1;
  const size_t beam_nodes = EulerBeam::mesh.get_nodes();

  ELFF_ASSERT(s.ints.size() == 5,
              "IBEulerBeamHuang::unpack_state(): invalid integer metadata.\n");
  ELFF_ASSERT(s.ints[0] == state_version,
              "IBEulerBeamHuang::unpack_state(): unsupported state version.\n");
  ELFF_ASSERT(static_cast<size_t>(s.ints[1]) == beam_nodes,
              "IBEulerBeamHuang::unpack_state(): node count mismatch.\n");

  const size_t expected_reals = 2 + 12 * beam_nodes +
                                static_cast<size_t>(last_T.size());
  ELFF_ASSERT(s.reals.size() == expected_reals,
              "IBEulerBeamHuang::unpack_state(): invalid real buffer size.\n");

  size_t k = 0;
  t = s.reals[k++];
  last_dt = s.reals[k++];
  unpack_matrix(s, k, X_init);
  unpack_matrix(s, k, X_n);
  unpack_matrix(s, k, X_nm1);
  unpack_matrix(s, k, V_init);
  for (Eigen::Index i = 0; i < last_T.size(); ++i) {
    last_T(i) = s.reals[k++];
  }

  time_iter = static_cast<size_t>(s.ints[2]);
  initial_velocity_pending = (s.ints[3] != 0);
  implicit_bending = (s.ints[4] != 0);

  EulerBeamInextensibleHuang::update_mesh();
  EBMeshToIBMeshCurrent();
  EBMeshToIBMeshNext();
}

} // namespace Models
} // namespace ELFF
