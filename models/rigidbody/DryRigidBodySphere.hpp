#pragma once

#include "elff/models/rigidbody/RigidContactBody.hpp"
#include "elff/models/runtime/Model.hpp"
#include "elff/models/rigidbody/RigidBodySphere.hpp"

namespace ELFF::Models {

/** @brief Runtime-owned dry sphere with analytic collision geometry.
 * Surface markers serve output and traction integration only.
 */
class DryRigidBodySphere : public Model, public RigidBodySphere,
                           public RigidContactBody {
public:
  using Vec3 = RigidBodyModel::Vec3;
  using Quat = std::array<real_t, 4>;

  DryRigidBodySphere(real_t radius, real_t density, size_t point_count,
      const Vec3& center = {}, const Vec3& velocity = {},
      const Quat& orientation = {1, 0, 0, 0}, const Vec3& angular_velocity = {});

  RigidBody3DModel& rigid_body() override { return *this; }
  const RigidBody3DModel& rigid_body() const override { return *this; }

  void pack_state(ModelState& state) const override;
  void unpack_state(const ModelState& state) override;
  void append_vtk_polydata(std::vector<IO::CXX::vtkPolyData>& datasets) const override;
};

} // namespace ELFF::Models
