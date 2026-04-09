#include "elff/models/beam/IBEulerBeamADDM.hpp"

namespace ELFF {
namespace Models {

void
IBEulerBeamADDM::EBMeshToIBMeshNext()
{
  EulerBeamInextensibleADDM::update_mesh();
  const size_t beam_nodes = EulerBeam::mesh.get_nodes();

  auto& ib_points = IBForceCoupled::mesh_next.GetPoints();
  auto& ib_velocity = IBForceCoupled::mesh_next.GetVelocity();

  auto& eb_points = EulerBeam::mesh.get_centerline();
  auto& eb_velocity = EulerBeam::mesh.get_centerline_velocity();

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
IBEulerBeamADDM::EBMeshToIBMeshCurrent()
{
  EulerBeamInextensibleADDM::update_mesh();
  const size_t beam_nodes = EulerBeam::mesh.get_nodes();

  auto& ib_points = IBForceCoupled::mesh.GetPoints();
  auto& ib_velocity = IBForceCoupled::mesh.GetVelocity();

  auto& eb_points = EulerBeam::mesh.get_centerline();
  auto& eb_velocity = EulerBeam::mesh.get_centerline_velocity();

  for (size_t ni = 0; ni < beam_nodes; ++ni) {
    ib_points[ni].x = eb_points[ni][0];
    ib_points[ni].y = eb_points[ni][1];
    ib_points[ni].z = eb_points[ni][2];

    ib_velocity[ni].x = eb_velocity[ni][0];
    ib_velocity[ni].y = eb_velocity[ni][1];
    ib_velocity[ni].z = eb_velocity[ni][2];
  }
}

IBEulerBeamADDM::IBEulerBeamADDM(real_t length,
                                 real_t EI,
                                 real_t mu,
                                 size_t nodes,
                                 EulerBeamBCs bcs,
                                 real_t r_penalty)
  : EulerBeamInextensibleADDM(length, EI, mu, nodes, bcs, r_penalty)
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
IBEulerBeamADDM::apply_initial_condition(EulerBeamMesh& mesh)
{
  EulerBeamInextensibleADDM::apply_initial_condition(mesh);
  EBMeshToIBMeshCurrent();
  EBMeshToIBMeshNext();
}

void
IBEulerBeamADDM::ComputeNextPoints(std::vector<IBMesh::IBVertex> force,
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

  EulerBeamInextensibleADDM::solve(dt, load);
  EBMeshToIBMeshNext();
}

void
IBEulerBeamADDM::pack_state(IBModelState& s) const
{
  static constexpr int64_t state_version = 1;
  const size_t beam_nodes = EulerBeam::mesh.get_nodes();

  s.ints.clear();
  s.reals.clear();
  s.bytes.clear();

  s.ints.push_back(state_version);
  s.ints.push_back(static_cast<int64_t>(beam_nodes));
  s.ints.push_back(static_cast<int64_t>(dof));
  s.ints.push_back(static_cast<int64_t>(elements));
  s.ints.push_back(static_cast<int64_t>(time_iter));
  s.ints.push_back(static_cast<int64_t>(have_prev_uniform_load));
  s.ints.push_back(static_cast<int64_t>(have_prev_nodal_load));
  s.ints.push_back(static_cast<int64_t>(nodal_load_prev.size()));

  const size_t collocation_dof = beam_nodes + elements;
  s.reals.reserve(1 + 18 * dof + 6 * collocation_dof + 3 +
                  3 * nodal_load_prev.size());

  s.reals.push_back(t);

  auto pack_vec = [&](const VectorXd& v) {
    for (Eigen::Index i = 0; i < v.size(); ++i)
      s.reals.push_back(v(i));
  };

  pack_vec(x);
  pack_vec(y);
  pack_vec(z);
  pack_vec(lambda_x);
  pack_vec(lambda_y);
  pack_vec(lambda_z);
  pack_vec(p);
  pack_vec(q);
  pack_vec(r);
  pack_vec(x_prev);
  pack_vec(y_prev);
  pack_vec(z_prev);
  pack_vec(vx_prev);
  pack_vec(vy_prev);
  pack_vec(vz_prev);
  pack_vec(ax_prev);
  pack_vec(ay_prev);
  pack_vec(az_prev);

  s.reals.push_back(load_prev[0]);
  s.reals.push_back(load_prev[1]);
  s.reals.push_back(load_prev[2]);

  for (const auto& load : nodal_load_prev) {
    s.reals.push_back(load[0]);
    s.reals.push_back(load[1]);
    s.reals.push_back(load[2]);
  }
}

void
IBEulerBeamADDM::unpack_state(const IBModelState& s)
{
  static constexpr int64_t state_version = 1;
  const size_t beam_nodes = EulerBeam::mesh.get_nodes();

  ELFF_ASSERT(s.ints.size() == 8,
              "IBEulerBeamADDM::unpack_state(): invalid integer metadata.\n");
  ELFF_ASSERT(s.ints[0] == state_version,
              "IBEulerBeamADDM::unpack_state(): unsupported state version.\n");
  ELFF_ASSERT(static_cast<size_t>(s.ints[1]) == beam_nodes,
              "IBEulerBeamADDM::unpack_state(): node count mismatch.\n");
  ELFF_ASSERT(static_cast<size_t>(s.ints[2]) == dof,
              "IBEulerBeamADDM::unpack_state(): dof mismatch.\n");
  ELFF_ASSERT(static_cast<size_t>(s.ints[3]) == elements,
              "IBEulerBeamADDM::unpack_state(): element count mismatch.\n");

  const size_t collocation_dof = beam_nodes + elements;
  const size_t nodal_load_size = static_cast<size_t>(s.ints[7]);
  const size_t expected_reals = 1 + 18 * dof + 6 * collocation_dof + 3 +
                                3 * nodal_load_size;
  ELFF_ASSERT(
    s.reals.size() == expected_reals,
    "IBEulerBeamADDM::unpack_state(): invalid real buffer size.\n");

  size_t k = 0;
  t = s.reals[k++];

  auto unpack_vec = [&](VectorXd& v) {
    for (Eigen::Index i = 0; i < v.size(); ++i)
      v(i) = s.reals[k++];
  };

  unpack_vec(x);
  unpack_vec(y);
  unpack_vec(z);
  unpack_vec(lambda_x);
  unpack_vec(lambda_y);
  unpack_vec(lambda_z);
  unpack_vec(p);
  unpack_vec(q);
  unpack_vec(r);
  unpack_vec(x_prev);
  unpack_vec(y_prev);
  unpack_vec(z_prev);
  unpack_vec(vx_prev);
  unpack_vec(vy_prev);
  unpack_vec(vz_prev);
  unpack_vec(ax_prev);
  unpack_vec(ay_prev);
  unpack_vec(az_prev);

  load_prev[0] = s.reals[k++];
  load_prev[1] = s.reals[k++];
  load_prev[2] = s.reals[k++];

  nodal_load_prev.resize(nodal_load_size);
  for (size_t i = 0; i < nodal_load_size; ++i) {
    nodal_load_prev[i][0] = s.reals[k++];
    nodal_load_prev[i][1] = s.reals[k++];
    nodal_load_prev[i][2] = s.reals[k++];
  }

  time_iter = static_cast<size_t>(s.ints[4]);
  have_prev_uniform_load = (s.ints[5] != 0);
  have_prev_nodal_load = (s.ints[6] != 0);

  EulerBeamInextensibleADDM::update_mesh();
  EBMeshToIBMeshCurrent();
  EBMeshToIBMeshNext();
}

} // namespace Models
} // namespace ELFF
