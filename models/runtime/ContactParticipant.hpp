#pragma once
#include "elff/contact/ContactSystem.hpp"

namespace ELFF::Models {
/** @brief A mechanical model participating in contact, independent of its solver.
 * Shapes are current world geometry, with stable local IDs. The runtime stamps
 * their owner IDs. Loads are integrated mechanical forces, never densities.
 * Geometry sampling and force mapping must be work conjugate. Body-pair eligibility
 * is specified per shape; self-contact is currently excluded by the evaluator.
 */
class ContactParticipant {
public:
  virtual ~ContactParticipant() = default;
  virtual std::vector<Contact::ShapeRecord> contact_shapes() const = 0;
  virtual void clear_contact_loads() = 0;
  virtual void add_contact_force(size_t shape, real_t coordinate, const Contact::Vec3& point,
                                 const Contact::Vec3& force) = 0;
  /** @brief Contact-only resultant; reference point is defined by each model. */
  virtual Contact::Wrench contact_wrench() const = 0;
  /** @brief Advance once with assembled loads held fixed; add gravity internally. */
  virtual void advance_contact(real_t dt, const Contact::Vec3& gravity) = 0;
  virtual real_t contact_time() const = 0;
  virtual size_t contact_iteration() const = 0;
  /** @brief Align the outer clock after accepted substeps, without resetting state. */
  virtual void synchronize_contact_time(real_t time) = 0;
  virtual void validate_contact_participation() const = 0;
};

/** @brief Optional split-step capability. Used only if every participant supports it. */
class SplitContactParticipant : public ContactParticipant {
public:
  virtual void predict_contact(real_t dt, const Contact::Vec3& gravity) = 0;
  virtual void finish_contact(real_t dt, const Contact::Vec3& gravity) = 0;
};
} // namespace ELFF::Models
