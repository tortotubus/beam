#pragma once
#include "elff/models/rigidbody/RigidBodyFibre.hpp"
#include "elff/models/rigidbody/RigidContactBody.hpp"
#include "elff/models/runtime/Model.hpp"

namespace ELFF::Models {
/** @brief Dry cylinder mechanics with an analytic capsule collision envelope.
 * The capsule axis is the full marker centerline; caps extend by diameter/2.
 */
class DryRigidBodyFibre : public Model, public RigidBodyFibre, public RigidContactBody {
public:
  using Quat = std::array<real_t, 4>;
  DryRigidBodyFibre(real_t length, real_t diameter, size_t nodes, real_t linear_density,
                    const Vec3& center = {}, const Vec3& velocity = {},
                    const Quat& orientation = {1, 0, 0, 0}, const Vec3& angular_velocity = {});
  RigidBody3DModel& rigid_body() override { return *this; }
  const RigidBody3DModel& rigid_body() const override { return *this; }
  void pack_state(ModelState& state) const override;
  void unpack_state(const ModelState& state) override;
  void append_vtk_polydata(std::vector<IO::CXX::vtkPolyData>& datasets) const override;
};
} // namespace ELFF::Models
