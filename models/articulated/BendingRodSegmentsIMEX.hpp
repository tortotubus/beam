#pragma once

#include "elff/models/articulated/BendingRodSegments.hpp"

namespace ELFF {
namespace Models {

using BendingRodSegmentsIMEXConfig = BendingRodSegmentsParameters;

/** Second-order extrapolated IMEX rod with three linear solves per timestep.
 *
 * Combine a full first-order step and two half steps as 2*refined-full,
 * then project onto unit tangents and transverse angular velocities.
 * Each stage linearizes bending, damping and velocity-dependent inertial
 * forces, with mass and geometry evaluated at that stage's starting state.
 * Stage tangent updates use exact rotations. The extrapolation cancels the
 * leading first-order error; two half steps alone would remain first order.
 * GetNext(force, dt) advances one fixed dt interval without adaptivity.
 * External forces are held constant during the interval. Second-order
 * fluid coupling would also require a time-accurate approximation of loads.
 * Large rotations can still require smaller steps for stability/accuracy.
 */
class BendingRodSegmentsIMEX : public BendingRodSegments {
public:
  explicit BendingRodSegmentsIMEX(const BendingRodSegmentsIMEXConfig &config);

  // Integrated endpoint forces, as in rhs(); ComputeNextPoints accepts
  // force densities and performs the quadrature conversion in the base.
  Eigen::VectorXd imexStep(const Eigen::VectorXd &state, double dt,
                           const std::vector<Vec3> &nodal_forces = {}) const;

  void pack_state(IBModelState &s) const override;
  void unpack_state(const IBModelState &s) override;

protected:
  Eigen::VectorXd integrateState(const Eigen::VectorXd &state, double dt,
                                 const std::vector<Vec3> &loads) const override;

private:
  Eigen::VectorXd firstOrderStep(const Eigen::VectorXd &state, double dt,
                                 const std::vector<Vec3> &nodal_forces) const;

  // Configuration-only translational mass coefficients. Each entry couples
  // the angular velocities of two segments through their center velocities.
  Eigen::MatrixXd translation_mass_;
};

} // namespace Models
} // namespace ELFF
