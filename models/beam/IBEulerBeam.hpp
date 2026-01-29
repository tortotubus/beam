#pragma once

#include "../ibm/immersedboundary.hpp"
#include "EulerBeamDynamicInextensibleMoM.hpp"

namespace ELFF {

/**
 * @brief The immersed boundary force-coupling class for @ref EulerBeam
 */
class IBEulerBeam
  : public IBForceCoupledStructureModel
  , public EulerBeamDynamicInextensibleMoM

{
private:
  void EBMeshToIBMeshNext()
  {
    // Update the E-B Mesh
    EulerBeamDynamicInextensibleMoM::update_mesh();

    // Get refernce to IB Mesh
    auto& ib_points = IBForceCoupledStructureModel::mesh_next.GetPoints();
    auto& ib_velocity = IBForceCoupledStructureModel::mesh_next.GetVelocity();

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
    auto& ib_points = IBForceCoupledStructureModel::mesh.GetPoints();
    auto& ib_velocity = IBForceCoupledStructureModel::mesh.GetVelocity();

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
    , IBForceCoupledStructureModel(nodes)
  {
    EBMeshToIBMeshCurrent();
    EBMeshToIBMeshNext();
  };

  void apply_initial_condition(EulerBeamMesh &mesh) override {
    EulerBeamDynamicInextensibleMoM::apply_initial_condition(mesh);
    EBMeshToIBMeshCurrent();
  }

  void ComputeNextPoints(std::vector<IBStructureMesh::IBVertex> force,
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

    auto& ib_points = IBForceCoupledStructureModel::mesh_next.GetPoints();
    auto& ib_velocity = IBForceCoupledStructureModel::mesh_next.GetVelocity();

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
};

}