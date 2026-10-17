#include "elff/models/articulated/BendingRodSegmentsIMEX.hpp"

#include <Eigen/Geometry>
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace ELFF {
namespace Models {
namespace {
using Basis = Eigen::Matrix<double, 3, 2>;

BendingRodSegmentsConfig sharedConfig(const BendingRodSegmentsParameters &config)
{
  BendingRodSegmentsConfig shared;
  static_cast<BendingRodSegmentsParameters &>(shared) = config;
  return shared;
}

Basis tangentBasis(const Vec3 &tangent)
{
  Eigen::Index axis;
  tangent.cwiseAbs().minCoeff(&axis);
  Basis basis;
  basis.col(0) = tangent.cross(Vec3::Unit(axis)).normalized();
  basis.col(1) = tangent.cross(basis.col(0));
  return basis;
}
} // namespace

BendingRodSegmentsIMEX::BendingRodSegmentsIMEX(
    const BendingRodSegmentsIMEXConfig &config)
  : BendingRodSegments(sharedConfig(config))
{
  Eigen::MatrixXd weights = Eigen::MatrixXd::Zero(n_, n_);
  for (int body = 0; body < n_; ++body) {
    for (int link = 0; link < n_; ++link) {
      const double lever = link < body ? h_ : (link == body ? h_ / 2 : 0.0);
      weights(body, link) = lever - (isFree() ? com_weights_[link] : 0.0);
    }
  }
  const Eigen::Map<const Eigen::VectorXd> masses(masses_.data(), n_);
  translation_mass_ = weights.transpose() * masses.asDiagonal() * weights;
}

Eigen::VectorXd BendingRodSegmentsIMEX::integrateState(
    const Eigen::VectorXd &state, double dt, const std::vector<Vec3> &loads) const
{
  return imexStep(state, dt, loads);
}

Eigen::VectorXd BendingRodSegmentsIMEX::imexStep(
    const Eigen::VectorXd &state, double dt,
    const std::vector<Vec3> &nodal_forces) const
{
  if (!std::isfinite(dt) || dt <= 0.0)
    throw std::invalid_argument("IMEX timestep must be finite and positive");
  // Richardson extrapolation cancels the O(dt^2) local error of the
  // first-order stage. Recompute geometry and forces at the half-step state.
  const Eigen::VectorXd full = firstOrderStep(state, dt, nodal_forces);
  const Eigen::VectorXd half = firstOrderStep(state, dt / 2, nodal_forces);
  const Eigen::VectorXd refined = firstOrderStep(half, dt / 2, nodal_forces);
  Eigen::VectorXd result = 2.0 * refined - full;
  if (!result.allFinite())
    throw std::runtime_error("non-finite extrapolated IMEX state");
  // Smooth projection preserves second order: the extrapolated state's
  // constraint defect is O(dt^3). Restore unit tangents and transverse spin.
  normalizeState(result);
  return result;
}

Eigen::VectorXd BendingRodSegmentsIMEX::firstOrderStep(
    const Eigen::VectorXd &state, double dt,
    const std::vector<Vec3> &nodal_forces) const
{
  const auto parts = unpack(state);
  // The existing articulated solve supplies the old-state acceleration, with
  // current bending and damping included. No duplicated force/constraint code.
  const Eigen::VectorXd derivative = rhs(0.0, state, nodal_forces);
  const int dofs = 2 * n_;
  std::vector<Basis> basis(n_), tangent_rate(n_);
  Eigen::VectorXd velocity(dofs), acceleration(dofs);
  for (int i = 0; i < n_; ++i) {
    basis[i] = tangentBasis(parts.tangents[i]);
    for (int axis = 0; axis < 2; ++axis)
      tangent_rate[i].col(axis) = basis[i].col(axis).cross(parts.tangents[i]);
    velocity.segment<2>(2 * i) = basis[i].transpose() * parts.omega[i];
    acceleration.segment<2>(2 * i) = basis[i].transpose() *
                                    derivative.segment<3>(omegaOffset() + 3 * i);
  }

  Eigen::MatrixXd mass(dofs, dofs);
  for (int i = 0; i < n_; ++i) {
    for (int j = 0; j < n_; ++j)
      mass.block<2, 2>(2 * i, 2 * j) = translation_mass_(i, j) *
                                      tangent_rate[i].transpose() * tangent_rate[j];
    mass.block<2, 2>(2 * i, 2 * i).diagonal().array() += inertias_[i];
  }

  Eigen::MatrixXd stiffness = Eigen::MatrixXd::Zero(dofs, dofs);
  Eigen::MatrixXd damping = Eigen::MatrixXd::Zero(dofs, dofs);
  for (int i = 1; i < n_; ++i) {
    const Vec3 &a = parts.tangents[i - 1];
    const Vec3 &b = parts.tangents[i];
    const Vec3 cross = a.cross(b);
    const double sine = cross.norm();
    const double cosine = std::clamp(a.dot(b), -1.0, 1.0);
    const double angle = std::atan2(sine, cosine);
    if (cosine < 0 && sine < 1e-8)
      throw std::runtime_error("IMEX bending linearization is undefined at antiparallel tangents");
    const double ratio = sine > 1e-6 ? angle / sine : 1.0 + sine * sine / 6.0;
    Vec3 mean = a + b;
    mean.normalize();
    const Mat3 projector = Mat3::Identity() - mean * mean.transpose();

    // Differentiate child torque -k * atan2(|a x b|,a.b) * (a x b)/|a x b|
    // with respect to each endpoint's two admissible rotation increments.
    for (int source : {i - 1, i}) {
      const double source_sign = source == i ? 1.0 : -1.0;
      for (int axis = 0; axis < 2; ++axis) {
        Vec3 da = Vec3::Zero(), db = Vec3::Zero();
        if (source == i - 1) da = tangent_rate[source].col(axis);
        else db = tangent_rate[source].col(axis);
        const Vec3 dcross = da.cross(b) + a.cross(db);
        const double dcosine = da.dot(b) + a.dot(db);
        double dratio;
        if (sine > 1e-6) {
          const double dsine = cross.dot(dcross) / sine;
          const double dangle = cosine * dsine - sine * dcosine;
          dratio = (dangle * sine - angle * dsine) / (sine * sine);
        } else {
          dratio = cross.dot(dcross) / 3.0;
        }
        const Vec3 dtorque = -bending_stiffness_[i - 1] *
                             (ratio * dcross + dratio * cross);
        const Vec3 damp = bending_damping_[i - 1] * source_sign *
                          projector * basis[source].col(axis);
        stiffness.block<2, 1>(2 * i, 2 * source + axis) -= basis[i].transpose() * dtorque;
        stiffness.block<2, 1>(2 * (i - 1), 2 * source + axis) += basis[i - 1].transpose() * dtorque;
        damping.block<2, 1>(2 * i, 2 * source + axis) += basis[i].transpose() * damp;
        damping.block<2, 1>(2 * (i - 1), 2 * source + axis) -= basis[i - 1].transpose() * damp;
      }
    }
  }

  // Linearize the velocity-dependent centrifugal force as well. With
  // omega perpendicular to t, omega x (omega x t) = -t |omega|^2.
  // Its generalized-force Jacobian is G_ij = 2 W_ij T_i^T t_j v_j^T.
  // Freezing this force explicitly can destabilize large-amplitude swings
  // even when the elastic modes are implicit.
  Eigen::MatrixXd inertial_jacobian(dofs, dofs);
  for (int i = 0; i < n_; ++i)
    for (int j = 0; j < n_; ++j)
      inertial_jacobian.block<2, 2>(2 * i, 2 * j) =
          2.0 * translation_mass_(i, j) *
          (tangent_rate[i].transpose() * parts.tangents[j]) *
          velocity.segment<2>(2 * j).transpose();

  // M (v_new-v)/dt = M a_old - dt K v_new + (G-D)(v_new-v).
  const Eigen::MatrixXd implicit_velocity = damping - inertial_jacobian;
  const Eigen::MatrixXd system = mass + dt * implicit_velocity + dt * dt * stiffness;
  const Eigen::VectorXd right = mass * (velocity + dt * acceleration) +
                                dt * implicit_velocity * velocity;
  const Eigen::VectorXd next_velocity = system.partialPivLu().solve(right);
  if (!next_velocity.allFinite() ||
      (system * next_velocity - right).norm() > 1e-9 * (1.0 + right.norm()))
    throw std::runtime_error("IMEX bending solve failed");

  Eigen::VectorXd result = state;
  for (int i = 0; i < n_; ++i) {
    const Vec3 omega = basis[i] * next_velocity.segment<2>(2 * i);
    const double speed = omega.norm();
    Vec3 tangent = parts.tangents[i];
    if (speed > 0)
      tangent = Eigen::AngleAxisd(dt * speed, omega / speed) * tangent;
    result.segment<3>(tangentOffset() + 3 * i) = tangent;
    result.segment<3>(omegaOffset() + 3 * i) = omega;
  }
  if (isFree()) {
    const Vec3 com_acceleration = derivative.segment<3>(velocityOffset());
    result.head<3>() = parts.center + dt * parts.center_velocity +
                       0.5 * dt * dt * com_acceleration;
    result.segment<3>(velocityOffset()) = parts.center_velocity + dt * com_acceleration;
  }
  if (!result.allFinite())
    throw std::runtime_error("non-finite state during IMEX integration");
  normalizeState(result);
  return result;
}

void BendingRodSegmentsIMEX::pack_state(IBModelState &s) const
{
  BendingRodSegments::pack_state(s);
  s.ints.push_back(2); // Second-order scheme tag; reject first-order/RK45 states.
}

void BendingRodSegmentsIMEX::unpack_state(const IBModelState &s)
{
  if (s.ints.size() != 5 || s.ints.back() != 2)
    throw std::invalid_argument("incompatible IMEX checkpoint");
  IBModelState common = s;
  common.ints.pop_back();
  BendingRodSegments::unpack_state(common);
}

} // namespace Models
} // namespace ELFF
