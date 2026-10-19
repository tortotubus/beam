#pragma once

#include "elff/config/config.hpp"
#include "elff/general/error.hpp"
#include <array>
#include <cstddef>
#include <vector>

namespace ELFF::Models {
/** Common nodal geometry and velocity, independent of the rod discretization.
 * Solvers publish their state here after initialization, solve and restore.
 * The material coordinate is uniformly spaced on [0, length].
 */
class BendingRodMesh {
public:
  using Vector = std::array<real_t, 3>;
  BendingRodMesh(size_t nodes, real_t length);
  size_t get_nodes() const { return s_.size(); }
  real_t get_length() const { return length_; }
  real_t get_ds() const { return ds_; }
  const std::vector<real_t>& get_curvilinear_axis() const { return s_; }
  real_t get_curvilinear_axis(size_t i) const {
    ELFF_ASSERT(i < s_.size(), "Rod mesh node index out of range.");
    return s_[i];
  }
  std::vector<Vector>& get_centerline() { return centerline_; }
  const std::vector<Vector>& get_centerline() const { return centerline_; }
  Vector& get_centerline(size_t i) {
    ELFF_ASSERT(i < centerline_.size(), "Rod mesh node index out of range.");
    return centerline_[i];
  }
  const Vector& get_centerline(size_t i) const {
    ELFF_ASSERT(i < centerline_.size(), "Rod mesh node index out of range.");
    return centerline_[i];
  }
  std::vector<Vector>& get_centerline_velocity() { return velocity_; }
  const std::vector<Vector>& get_centerline_velocity() const { return velocity_; }
  Vector& get_centerline_velocity(size_t i) {
    ELFF_ASSERT(i < velocity_.size(), "Rod mesh node index out of range.");
    return velocity_[i];
  }
  const Vector& get_centerline_velocity(size_t i) const {
    ELFF_ASSERT(i < velocity_.size(), "Rod mesh node index out of range.");
    return velocity_[i];
  }

protected:
  real_t length_, ds_;
  std::vector<real_t> s_;
  std::vector<Vector> centerline_, velocity_;
};
} // namespace ELFF::Models
