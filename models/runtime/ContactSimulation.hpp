#pragma once

#include "elff/contact/ContactSystem.hpp"
#include "elff/models/runtime/ModelRuntime.hpp"
#include "elff/models/runtime/ContactParticipant.hpp"
#include "elff/contact/CollisionGeometry.hpp"

namespace ELFF::Models {

/** @brief Runtime-owned serial contact simulation, supporting rigid and deformable participants.
 * Owns registration, fixed planes, contact evaluation and synchronized subcycling.
 * Rigid-only runs use split stepping; whole-step participants use held loads.
 * Models remain runtime-owned. Configure before sealing the runtime;
 * advancement failures do not roll back. Reconstruct settings for restart.
 */
class ContactSimulation : public RuntimeComponent {
public:
  using RuntimeComponent::RuntimeComponent;
  /** @brief Configure contact search, normal law, substep limit and gravity.
   * Configuration and registration are frozen by runtime sealing.
   */
  void configure(const Contact::ContactGrid& grid, Contact::ContactParameters parameters,
                 real_t max_substep, const Contact::Vec3& gravity = {});
  void register_body(ModelId id);
  /** @brief Add a fixed plane; return its stable zero-based insertion ID. */
  size_t add_plane(const Contact::Vec3& point, const Contact::Vec3& normal);
  /** @brief Resolve local participants and validate synchronized clocks; does not advance. */
  int seal() override;
  /** @brief Evaluate loads without advancement; requires a sealed runtime. */
  void evaluate();
  /** @brief Advance a finite nonnegative interval; requires a sealed runtime. */
  void advance(real_t interval);
  real_t time() const;
  real_t max_substep() const { return max_substep_; }
  bool uses_split_stepping() const { return split_stepping_; }
  size_t last_substeps() const { return last_substeps_; }

  const std::vector<ModelId>& model_ids() const { return ids_; }
  /** @brief Contact-only resultant at the last force evaluation; reference is model-defined. */
  Contact::Wrench wrench(ModelId id) const;
  const std::vector<Contact::ContactRecord>& contacts() const;
  const std::vector<Contact::Wrench>& plane_reactions() const;

private:
  void require_mutable() const;
  void require_sealed() const;
  void verify_synchronized() const;
  std::vector<ModelId> ids_;
  std::vector<ContactParticipant*> bodies_;
  std::vector<SplitContactParticipant*> split_bodies_;
  bool split_stepping_ = false;
  std::vector<Contact::PlaneGeometry> planes_;
  std::unique_ptr<Contact::ContactSystem> system_;
  real_t max_substep_ = 0;
  Contact::Vec3 gravity_{};
  size_t last_substeps_ = 0;
  bool sealed_ = false;
};

} // namespace ELFF::Models
