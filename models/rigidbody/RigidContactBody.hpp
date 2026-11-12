#pragma once
#include "elff/models/rigidbody/RigidBody3DModel.hpp"
#include "elff/models/runtime/ContactParticipant.hpp"

namespace ELFF::Models {
/** @brief Rigid implementation of the generic contact participant capability.
 * IBM adapters do not opt in until coupling-mesh synchronization is implemented.
 */
class RigidContactBody : public SplitContactParticipant {
public:
  virtual RigidBody3DModel& rigid_body() = 0;
  virtual const RigidBody3DModel& rigid_body() const = 0;
  std::vector<Contact::ShapeRecord> contact_shapes() const override;
  void clear_contact_loads() override { load_ = {}; }
  void add_contact_force(size_t shape, real_t coordinate, const Contact::Vec3& point,
                         const Contact::Vec3& force) override;
  Contact::Wrench contact_wrench() const override { return load_; }
  void advance_contact(real_t dt, const Contact::Vec3& gravity) override;
  void predict_contact(real_t dt, const Contact::Vec3& gravity) override;
  void finish_contact(real_t dt, const Contact::Vec3& gravity) override;
  real_t contact_time() const override { return rigid_body().time(); }
  size_t contact_iteration() const override { return rigid_body().time_iteration(); }
  void synchronize_contact_time(real_t time) override;
  void validate_contact_participation() const override;

private:
  Contact::Wrench load_;
};
} // namespace ELFF::Models
