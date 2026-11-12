#pragma once
#include "elff/models/articulated/BendingRodImplicit.hpp"
#include "elff/models/runtime/ContactParticipant.hpp"

namespace ELFF::Models {
/** @brief Contact adapter for an implicit articulated rod.
 * One world capsule per rigid link. Contact loads are work-conjugate endpoint
 * forces. Radius controls only collision geometry; it does not change inertia.
 * Body-pair contact is opt-in; self-contact and IBM coupling remain disabled.
 */
class DryBendingRodImplicit : public Model, public BendingRodImplicit, public ContactParticipant {
public:
  DryBendingRodImplicit(const BendingRodImplicitConfig& config, real_t radius,
                        bool body_contacts = false);
  using BendingRodImplicit::solve;
  void solve(real_t dt, const std::vector<Vector>& density) override;
  void solve_integrated(real_t dt, const std::vector<Vector>& forces) override;
  void apply_initial_condition(const BendingRodMesh& mesh) override;
  std::vector<Contact::ShapeRecord> contact_shapes() const override;
  void clear_contact_loads() override;
  void add_contact_force(size_t shape, real_t coordinate, const Contact::Vec3& point,
                         const Contact::Vec3& force) override;
  /** @brief Contact-only resultant, moment about the current first mesh node. */
  Contact::Wrench contact_wrench() const override;
  const std::vector<Vector>& contact_nodal_forces() const { return loads_; }
  real_t collision_radius() const { return radius_; }
  bool body_contacts_enabled() const { return body_contacts_; }
  void advance_contact(real_t dt, const Contact::Vec3& gravity) override;
  real_t contact_time() const override { return time_; }
  size_t contact_iteration() const override { return iteration_; }
  void synchronize_contact_time(real_t time) override;
  void validate_contact_participation() const override;
  void pack_state(ModelState& state) const override;
  void unpack_state(const ModelState& state) override;
  void append_vtk_polydata(std::vector<IO::CXX::vtkPolyData>& datasets) const override;

private:
  void validate_step(real_t dt) const;
  real_t radius_;
  bool body_contacts_;
  real_t time_ = 0;
  size_t iteration_ = 0;
  std::vector<Vector> loads_;
};
} // namespace ELFF::Models
