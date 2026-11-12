#pragma once

#include "elff/contact/BroadPhaseLinkedCell.hpp"
#include "elff/contact/CollisionGeometry.hpp"
#include <optional>

namespace ELFF::Contact {

struct ContactParameters {
  real_t normal_stiffness;
  real_t normal_damping;
};

struct ContactGrid {
  Vec3 origin;
  real_t cell_width;
  std::array<size_t, 3> dimensions;
};

/** @brief One world shape owned by a model, with rigid or endpoint kinematics.
 * Identity is the pair (id, local_shape). Shape indices are not model indices.
 */
struct ShapeRecord {
  ShapeId id;
  CollisionShape geometry;
  Vec3 center;
  Vec3 velocity;
  Vec3 angular_velocity;
  size_t local_shape = 0;
  /** @brief Optional linearly interpolated centerline velocities for deformable capsules. */
  std::optional<std::array<Vec3, 2>> endpoint_velocities;
  /** @brief Optional persistent material-node identities for shared-endpoint deduplication. */
  std::optional<std::array<size_t, 2>> endpoint_features;
  bool body_contacts = true;
};

/** @brief One overlapping pair, including contacts whose force is clamped to zero.
 * For a plane contact, b is its insertion index. Otherwise b is a body ID.
 */
struct ContactRecord {
  ShapeId a;
  ShapeId b;
  bool plane;
  ContactPoint geometry;
  Vec3 force_on_a;
  size_t shape_a = 0;
  size_t shape_b = 0;
};

/** @brief Frictionless shape contact evaluation; never advances structural state.
 * Each evaluation replaces diagnostics. Shape wrenches follow input order;
 * plane reactions follow plane order, with moments about each plane's point.
 * No history or persistent geometry snapshots are maintained.
 */
class ContactSystem {
public:
  ContactSystem(const ContactGrid& grid, ContactParameters parameters);
  void evaluate(const std::vector<ShapeRecord>& bodies,
                const std::vector<PlaneGeometry>& planes);
  const std::vector<Wrench>& shape_wrenches() const { return shape_wrenches_; }
  const std::vector<Wrench>& plane_reactions() const { return plane_reactions_; }
  const std::vector<ContactRecord>& contacts() const { return contacts_; }

private:
  ContactParameters parameters_;
  BroadPhaseLinkedCell broad_phase_;
  std::vector<Wrench> shape_wrenches_;
  std::vector<Wrench> plane_reactions_;
  std::vector<ContactRecord> contacts_;
};

} // namespace ELFF::Contact
