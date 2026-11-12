#include "elff/models/runtime/ContactSimulation.hpp"
#include <algorithm>
#include <cmath>
#include <limits>

namespace ELFF::Models {

void ContactSimulation::require_mutable() const {
  ELFF_VERIFY(runtime_.mutable_registration() && !sealed_ && runtime_.ranks() == 1,
              "Contact registration requires an unsealed serial runtime.");
}

void ContactSimulation::configure(const Contact::ContactGrid& grid,
                                  Contact::ContactParameters parameters, real_t max_substep,
                                  const Contact::Vec3& gravity) {
  require_mutable();
  ELFF_VERIFY(std::isfinite(max_substep) && max_substep > 0,
              "Maximum substep must be finite and positive.");
  for (auto g : gravity)
    ELFF_VERIFY(std::isfinite(g), "Gravity must be finite.");
  auto system = std::make_unique<Contact::ContactSystem>(grid, parameters);
  system_ = std::move(system);
  max_substep_ = max_substep;
  gravity_ = gravity;
}

void ContactSimulation::register_body(ModelId id) {
  require_mutable();
  const auto* entry = runtime_.entry(id);
  ELFF_VERIFY(entry && entry->owner_pid == runtime_.rank() &&
                  dynamic_cast<ContactParticipant*>(entry->model.get()) &&
                  std::find(ids_.begin(), ids_.end(), id) == ids_.end(),
              "Contact registration requires a unique, local mechanical model ID.");
  ids_.push_back(id);
}

size_t ContactSimulation::add_plane(const Contact::Vec3& point, const Contact::Vec3& normal) {
  require_mutable();
  planes_.push_back(Contact::make_plane(point, normal));
  return planes_.size() - 1;
}

int ContactSimulation::seal() {
  if (sealed_)
    return 0;
  ELFF_VERIFY(system_ && runtime_.ranks() == 1,
              "Configure contact before sealing a serial runtime.");
  std::vector<ContactParticipant*> bodies;
  for (auto id : ids_) {
    const auto* entry = runtime_.entry(id);
    auto* participant = dynamic_cast<ContactParticipant*>(runtime_.model(id));
    auto* body = participant;
    ELFF_VERIFY(entry && entry->owner_pid == runtime_.rank() && body,
                "Contact body must remain a local mechanical model.");
    body->validate_contact_participation();
    bodies.push_back(body);
  }
  bodies_ = std::move(bodies);
  split_bodies_.clear();
  split_stepping_ = true;
  for (auto* body : bodies_) {
    auto* split = dynamic_cast<SplitContactParticipant*>(body);
    split_bodies_.push_back(split);
    split_stepping_ = split_stepping_ && split;
  }
  verify_synchronized();
  sealed_ = true;
  return 0;
}

void ContactSimulation::evaluate() {
  require_sealed();
  std::vector<Contact::ShapeRecord> records;
  for (size_t i = 0; i < bodies_.size(); ++i) {
    auto shapes = bodies_[i]->contact_shapes();
    ELFF_VERIFY(!shapes.empty(), "Contact participant must expose at least one shape.");
    for (auto& shape : shapes) {
      shape.id = static_cast<Contact::ShapeId>(ids_[i]);
      records.push_back(shape);
    }
    bodies_[i]->clear_contact_loads();
  }
  system_->evaluate(records, planes_);
  auto participant = [&](Contact::ShapeId id) {
    const auto found = std::find(ids_.begin(), ids_.end(), static_cast<ModelId>(id));
    ELFF_VERIFY(found != ids_.end(), "Unknown contact owner.");
    return bodies_[static_cast<size_t>(found - ids_.begin())];
  };
  for (const auto& contact : system_->contacts()) {
    participant(contact.a)->add_contact_force(contact.shape_a, contact.geometry.parameter_a,
                                              contact.geometry.point, contact.force_on_a);
    if (!contact.plane) {
      auto opposite = contact.force_on_a;
      for (auto& v : opposite)
        v = -v;
      participant(contact.b)->add_contact_force(contact.shape_b, contact.geometry.parameter_b,
                                                contact.geometry.point, opposite);
    }
  }
}

Contact::Wrench ContactSimulation::wrench(ModelId id) const {
  require_sealed();
  const auto found = std::find(ids_.begin(), ids_.end(), id);
  ELFF_VERIFY(found != ids_.end(), "Unregistered contact participant.");
  return bodies_[static_cast<size_t>(found - ids_.begin())]->contact_wrench();
}

const std::vector<Contact::ContactRecord>& ContactSimulation::contacts() const {
  ELFF_VERIFY(system_, "Contact is not configured.");
  return system_->contacts();
}

const std::vector<Contact::Wrench>& ContactSimulation::plane_reactions() const {
  ELFF_VERIFY(system_, "Contact is not configured.");
  return system_->plane_reactions();
}

void ContactSimulation::require_sealed() const {
  ELFF_VERIFY(runtime_.sealed() && sealed_,
              "Seal runtime before contact evaluation or advancement.");
}

void ContactSimulation::verify_synchronized() const {
  ELFF_VERIFY(!bodies_.empty(), "Contact simulation requires at least one registered body.");
  const auto* first = bodies_.front();
  ELFF_VERIFY(std::isfinite(first->contact_time()), "Participant clock must be finite.");
  for (const auto* body : bodies_)
    ELFF_VERIFY(body->contact_time() == first->contact_time() &&
                    body->contact_iteration() == first->contact_iteration(),
                "Participant clocks must be synchronized.");
}

real_t ContactSimulation::time() const {
  require_sealed();
  verify_synchronized();
  return bodies_.front()->contact_time();
}

void ContactSimulation::advance(real_t interval) {
  ELFF_VERIFY(std::isfinite(interval) && interval >= 0, "Interval must be finite and nonnegative.");
  require_sealed();
  last_substeps_ = 0;
  if (interval == 0)
    return;
  verify_synchronized();
  const real_t start = time();
  ELFF_VERIFY(std::isfinite(start + interval) && start + interval > start,
              "Interval must advance finite time.");
  const real_t count = std::max(real_t{1}, std::ceil(interval / max_substep_));
  // Strict upper bound avoids undefined conversion when size_t max rounds up.
  ELFF_VERIFY(std::isfinite(count) && count >= 1 &&
                  count < static_cast<real_t>(std::numeric_limits<size_t>::max()),
              "Substep count overflow.");
  size_t steps = static_cast<size_t>(count);
  real_t h = interval / static_cast<real_t>(steps);
  if (h > max_substep_) {
    ELFF_VERIFY(steps < std::numeric_limits<size_t>::max() - 1, "Substep count overflow.");
    h = interval / static_cast<real_t>(++steps);
  }
  ELFF_VERIFY(h > 0 && std::isfinite(h) && h <= max_substep_ && h / 2 > 0,
              "Substep is unrepresentable.");
  for (const auto* body : bodies_)
    ELFF_VERIFY(steps <= std::numeric_limits<size_t>::max() - body->contact_iteration(),
                "Participant iteration overflow.");
  // Whole-step participants retain the loads evaluated at the previous step's
  // end. Always refresh on entry: callers may have changed or restored state.
  if (!split_stepping_)
    evaluate();
  for (size_t m = 0; m < steps; ++m) {
    if (split_stepping_) {
      // The final kick changes velocities, so damping loads must be refreshed.
      evaluate();
      for (auto* body : split_bodies_)
        body->predict_contact(h, gravity_);
      evaluate();
      for (auto* body : split_bodies_)
        body->finish_contact(h, gravity_);
    } else {
      // Every model sees loads from the same configuration. Structural solvers
      // own their internal integration; contacts remain explicitly lagged.
      for (auto* body : bodies_)
        body->advance_contact(h, gravity_);
      evaluate();
    }
    ++last_substeps_;
  }
  // Avoid accumulated clock roundoff at an outer synchronization boundary.
  for (auto* body : bodies_)
    body->synchronize_contact_time(start + interval);
}

} // namespace ELFF::Models
