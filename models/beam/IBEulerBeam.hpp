#pragma once

#include "models/ibm/IBForceCoupled.hpp"
#include "models/beam/EulerBeamDynamicInextensibleMoM.hpp"


namespace ELFF {
namespace Models {

/**
 * @brief The immersed boundary force-coupling class for @ref EulerBeam
 */
class IBEulerBeam
  : public IBForceCoupled
  , public EulerBeamDynamicInextensibleMoM
{
private:
  void EBMeshToIBMeshNext()
  {
    // Update the E-B Mesh
    EulerBeamDynamicInextensibleMoM::update_mesh();

    // Get refernce to IB Mesh
    auto& ib_points = IBForceCoupled::mesh_next.GetPoints();
    auto& ib_velocity = IBForceCoupled::mesh_next.GetVelocity();

    // Get reference to E-B Mesh
    auto& eb_points = EulerBeam::mesh.get_centerline();
    auto& eb_velocity = EulerBeam::mesh.get_centerline_velocity();

    // Copy
    for (size_t ni = 0; ni < nodes; ni++) {
      ib_points[ni].x = eb_points[ni][0];
      ib_points[ni].y = eb_points[ni][1];
      ib_points[ni].z = eb_points[ni][2];

      ib_velocity[ni].x = eb_velocity[ni][0];
      ib_velocity[ni].y = eb_velocity[ni][1];
      ib_velocity[ni].z = eb_velocity[ni][2];
    }
  }

  void EBMeshToIBMeshCurrent()
  {
    // Update the E-B Mesh
    EulerBeamDynamicInextensibleMoM::update_mesh();

    // Get refernce to IB Mesh
    auto& ib_points = IBForceCoupled::mesh.GetPoints();
    auto& ib_velocity = IBForceCoupled::mesh.GetVelocity();

    // Get reference to E-B Mesh
    auto& eb_points = EulerBeam::mesh.get_centerline();
    auto& eb_velocity = EulerBeam::mesh.get_centerline_velocity();

    // Copy
    for (size_t ni = 0; ni < nodes; ni++) {
      ib_points[ni].x = eb_points[ni][0];
      ib_points[ni].y = eb_points[ni][1];
      ib_points[ni].z = eb_points[ni][2];

      ib_velocity[ni].x = eb_velocity[ni][0];
      ib_velocity[ni].y = eb_velocity[ni][1];
      ib_velocity[ni].z = eb_velocity[ni][2];
    }
  }

public:
  IBEulerBeam(real_t length,
              real_t EI,
              real_t mu,
              size_t nodes,
              EulerBeamBCs bcs,
              real_t r_penalty)
    : EulerBeamDynamicInextensibleMoM(length, EI, mu, nodes, bcs, r_penalty)
    , IBForceCoupled(nodes)
  {
    EBMeshToIBMeshCurrent();
    EBMeshToIBMeshNext();
  };

  void apply_initial_condition(EulerBeamMesh& mesh) override
  {
    EulerBeamDynamicInextensibleMoM::apply_initial_condition(mesh);
    EBMeshToIBMeshCurrent();
  }

  void ComputeNextPoints(std::vector<IBMesh::IBVertex> force,
                         real_t dt) override
  {
    ELFF_ASSERT(force.size() == nodes, "Force array must match node count.\n");

    std::vector<std::array<real_t, 3>> load(nodes);
    for (size_t ni = 0; ni < nodes; ni++) {
      load[ni][0] = force[ni].x;
      load[ni][1] = force[ni].y;
      load[ni][2] = force[ni].z;
    }

    EulerBeamDynamicInextensibleMoM::solve(dt, load);

    auto& ib_points = IBForceCoupled::mesh_next.GetPoints();
    auto& ib_velocity = IBForceCoupled::mesh_next.GetVelocity();

    auto& eb_points = EulerBeamDynamicInextensibleMoM::mesh.get_centerline();
    auto& eb_velocity =
      EulerBeamDynamicInextensibleMoM::mesh.get_centerline_velocity();

    for (size_t ni = 0; ni < nodes; ni++) {
      ib_points[ni].x = eb_points[ni][0];
      ib_points[ni].y = eb_points[ni][1];
      ib_points[ni].z = eb_points[ni][2];

      ib_velocity[ni].x = eb_velocity[ni][0];
      ib_velocity[ni].y = eb_velocity[ni][1];
      ib_velocity[ni].z = eb_velocity[ni][2];
    }
  }

  void pack_state(IBModelState &s) const override {
    static constexpr int64_t state_version = 1;

    s.ints.clear();
    s.reals.clear();
    s.bytes.clear();

    s.ints.push_back(state_version);
    s.ints.push_back(static_cast<int64_t>(nodes));
    s.ints.push_back(static_cast<int64_t>(ndof));
    s.ints.push_back(static_cast<int64_t>(ndof_l));
    s.ints.push_back(static_cast<int64_t>(time_iter));

    s.reals.reserve(1 + 3 * ndof + ndof_l);
    s.reals.push_back(t);

    auto pack_vec = [&](const VectorXd& v) {
      for (Eigen::Index i = 0; i < v.size(); ++i)
        s.reals.push_back(v(i));
    };

    pack_vec(u);
    pack_vec(lambda);
    pack_vec(v_prev);
    pack_vec(a_prev);
  }

  void unpack_state(const IBModelState& s) override
  {
    static constexpr int64_t state_version = 1;

    ELFF_ASSERT(s.ints.size() == 5,
                "IBEulerBeam::unpack_state(): invalid integer metadata.\n");
    ELFF_ASSERT(s.ints[0] == state_version,
                "IBEulerBeam::unpack_state(): unsupported state version.\n");
    ELFF_ASSERT(static_cast<size_t>(s.ints[1]) == nodes,
                "IBEulerBeam::unpack_state(): node count mismatch.\n");
    ELFF_ASSERT(static_cast<size_t>(s.ints[2]) == ndof,
                "IBEulerBeam::unpack_state(): ndof mismatch.\n");
    ELFF_ASSERT(static_cast<size_t>(s.ints[3]) == ndof_l,
                "IBEulerBeam::unpack_state(): constraint dof mismatch.\n");

    const size_t expected_reals = 1 + 3 * ndof + ndof_l;
    ELFF_ASSERT(s.reals.size() == expected_reals,
                "IBEulerBeam::unpack_state(): invalid real buffer size.\n");

    size_t k = 0;
    t = s.reals[k++];

    auto unpack_vec = [&](VectorXd& v) {
      for (Eigen::Index i = 0; i < v.size(); ++i)
        v(i) = s.reals[k++];
    };

    unpack_vec(u);
    unpack_vec(lambda);
    unpack_vec(v_prev);
    unpack_vec(a_prev);

    time_iter = static_cast<size_t>(s.ints[4]);

    // Reconstruct step history and all geometric/IB views from the restored
    // converged state.
    u_prev = u;
    u_prev_prev = u;
    EulerBeamDynamicInextensibleMoM::update_mesh();
    EBMeshToIBMeshCurrent();
    EBMeshToIBMeshNext();
  }
};

} // namespace Models 
} // namespace ELFF 
