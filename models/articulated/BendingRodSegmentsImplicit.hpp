#pragma once

#include "elff/models/articulated/BendingRodSegments.hpp"

namespace ELFF {
namespace Models {

struct BendingRodSegmentsImplicitConfig : BendingRodSegmentsParameters {
  double spectral_radius = 0.8; // [0,1]; smaller values damp fast linear modes.
  // Infinity norm of the mass-scaled momentum residual, in velocity units.
  double nonlinear_tolerance = 1e-9;
  int max_iterations = 20;
};

/** Fully implicit second-order generalized-alpha structural integrator.
 *
 * Each step uses exponential sphere coordinates
 * t(u)=cos(|u|)*t_old+sinc(|u|)*B*u, with two coordinates per segment.
 * Newmark kinematics and the complete nonlinear acceleration equation are
 * solved together in these coordinates. Unit tangents and transverse angular
 * velocities follow from the chart, without an endpoint projection.
 * Algorithmic accelerations are transported when recentering the charts.
 * There are 2*links unknowns (plus three COM unknowns for a free rod).
 * Newton uses a dense Jacobian assembled from local autodiff derivatives and backtracking.
 *
 * GetNext(force,dt) takes one fixed step; failure throws without advancing.
 * Loads supplied by the IB caller are held constant over that step: this does
 * not implement an implicit fluid solve. Large rotations/stiffness can require
 * smaller steps; linear generalized-alpha stability is not a nonlinear bound.
 * Energy need not decrease each step, even with numerical damping enabled.
 * Changing the timestep restarts acceleration history from the physical RHS.
 */
class BendingRodSegmentsImplicit : public BendingRodSegments {
public:
  explicit BendingRodSegmentsImplicit(const BendingRodSegmentsImplicitConfig &config);
  void setState(const Eigen::VectorXd &state) override;
  void pack_state(IBModelState &s) const override;
  void unpack_state(const IBModelState &s) override;

protected:
  Eigen::VectorXd integrateState(const Eigen::VectorXd &state, double dt,
                                const std::vector<Vec3> &loads) const override;

private:
  Eigen::MatrixXd tangent_mass_;
  double spectral_radius_, tolerance_;
  int max_iterations_;
  // Ambient tangent second derivatives, preceded by COM acceleration if free.
  mutable Eigen::VectorXd acceleration_;
  mutable double previous_dt_ = 0.0;
};

} // namespace Models
} // namespace ELFF
