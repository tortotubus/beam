#include "elff/models/articulated/BendingRodSegments.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace ELFF {
namespace Models {
namespace {

std::size_t checkedNodeCount(int links) {
  if (links < 1) {
    throw std::invalid_argument("links must be positive");
  }
  return static_cast<std::size_t>(links) + 1;
}

using Mat2 = Eigen::Matrix2d;
using Mat32 = Eigen::Matrix<double, 3, 2>;

Mat3 crossMatrix(const Vec3 &value) {
  Mat3 result;
  result << 0.0, -value.z(), value.y(), value.z(), 0.0, -value.x(), -value.y(),
      value.x(), 0.0;
  return result;
}

Mat32 tangentBasis(const Vec3 &tangent) {
  Eigen::Index least_aligned;
  tangent.cwiseAbs().minCoeff(&least_aligned);
  const Vec3 reference = Vec3::Unit(static_cast<int>(least_aligned));
  Mat32 basis;
  basis.col(0) = tangent.cross(reference).normalized();
  basis.col(1) = tangent.cross(basis.col(0));
  return basis;
}

template <typename Value>
std::vector<Value> configuredValues(const std::vector<Value> &supplied,
                                    int count, const Value &fallback,
                                    const char *name) {
  if (supplied.empty()) {
    return std::vector<Value>(static_cast<std::size_t>(count), fallback);
  }
  if (static_cast<int>(supplied.size()) != count) {
    throw std::invalid_argument(std::string(name) + " has the wrong size");
  }
  return supplied;
}

} // namespace

BendingRodSegments::BendingRodSegments(const BendingRodSegmentsConfig &config)
    : IBForceCoupled(checkedNodeCount(config.links)), n_(config.links), length_(config.length), total_mass_(config.total_mass),
      boundary_(config.boundary), gravity_(config.gravity), pin_(config.pin),
      integration_step_(config.integration_step), adaptive_(config.adaptive),
      relative_tolerance_(config.relative_tolerance),
      absolute_tolerance_(config.absolute_tolerance) {
  validateNonnegative(integration_step_, "integration_step");
  if (!std::isfinite(relative_tolerance_) || relative_tolerance_ <= 0.0 ||
      !std::isfinite(absolute_tolerance_) || absolute_tolerance_ <= 0.0) {
    throw std::invalid_argument("integration tolerances must be finite and positive");
  }
  if (n_ < 1 || !(length_ > 0.0) || !(total_mass_ > 0.0) ||
      !std::isfinite(length_) || !std::isfinite(total_mass_)) {
    throw std::invalid_argument(
        "links, length, and total_mass must be positive");
  }
  if (!std::isfinite(config.ei) || config.ei < 0.0 || !gravity_.allFinite() ||
      !pin_.allFinite()) {
    throw std::invalid_argument(
        "beam parameters must be finite and nonnegative");
  }
  h_ = length_ / static_cast<double>(n_);

  masses_ = configuredValues(config.masses, n_, total_mass_ / n_, "masses");
  total_mass_ = 0.0;
  for (double mass : masses_) {
    if (!std::isfinite(mass) || mass <= 0.0) {
      throw std::invalid_argument("masses must be finite and positive");
    }
    total_mass_ += mass;
  }

  inertias_.resize(static_cast<std::size_t>(n_));
  if (config.inertias.empty()) {
    for (int i = 0; i < n_; ++i) {
      inertias_[i] = masses_[i] * h_ * h_ / 12.0;
    }
  } else {
    inertias_ = configuredValues(config.inertias, n_, 0.0, "inertias");
  }
  for (double inertia : inertias_) {
    if (!std::isfinite(inertia) || inertia <= 0.0) {
      throw std::invalid_argument(
          "transverse inertias must be finite and positive");
    }
  }

  const int joints = n_ - 1;
  bending_stiffness_ = configuredValues(config.bending_stiffness, joints,
                                        config.ei / h_, "bending_stiffness");
  bending_damping_ =
      configuredValues(config.bending_damping, joints, 0.0, "bending_damping");
  for (int i = 0; i < joints; ++i) {
    validateNonnegative(bending_stiffness_[i], "bending_stiffness");
    validateNonnegative(bending_damping_[i], "bending_damping");
  }

  com_weights_.resize(n_);
  double downstream_mass = 0.0;
  for (int i = n_ - 1; i >= 0; --i) {
    downstream_mass += masses_[i];
    com_weights_[i] = h_ * (downstream_mass - 0.5 * masses_[i]) / total_mass_;
  }
  std::vector<real_t> measures(static_cast<std::size_t>(n_) + 1, h_);
  measures.front() *= 0.5;
  measures.back() *= 0.5;
  SetNodalMeasures(measures);
  setState(initialState());
}

int BendingRodSegments::links() const { return n_; }
int BendingRodSegments::stateSize() const { return 6 * n_ + (isFree() ? 6 : 0); }
bool BendingRodSegments::isFree() const { return boundary_ == Boundary::Free; }

void BendingRodSegments::setState(const Eigen::VectorXd &state) {
  unpack(state); // Validate before modifying the owned state.
  Eigen::VectorXd normalized = state;
  normalizeState(normalized);
  state_ = normalized;
  std::fill(mesh.GetForces().begin(), mesh.GetForces().end(), IBMesh::IBVertex{0, 0, 0});
  syncMeshes();
}

void BendingRodSegments::syncMeshes() {
  const auto parts = unpack(state_);
  const Eigen::MatrixXd positions = nodes(state_);
  std::vector<Vec3> velocities(static_cast<std::size_t>(n_) + 1, Vec3::Zero());
  Vec3 weighted_velocity = Vec3::Zero();
  for (int i = 0; i < n_; ++i) {
    velocities[i + 1] = velocities[i] + parts.omega[i].cross(h_ * parts.tangents[i]);
    weighted_velocity += 0.5 * masses_[i] * (velocities[i] + velocities[i + 1]);
  }
  const Vec3 shift = isFree() ? Vec3(parts.center_velocity - weighted_velocity / total_mass_)
                              : Vec3::Zero();
  for (int i = 0; i <= n_; ++i) {
    const Vec3 velocity = velocities[i] + shift;
    mesh.GetPoints()[i] = {static_cast<real_t>(positions(i, 0)),
                           static_cast<real_t>(positions(i, 1)),
                           static_cast<real_t>(positions(i, 2))};
    mesh.GetVelocity()[i] = {static_cast<real_t>(velocity.x()),
                             static_cast<real_t>(velocity.y()),
                             static_cast<real_t>(velocity.z())};
  }
  CopyCurrentToNext();
}

void BendingRodSegments::ComputeNextPoints(std::vector<IBMesh::IBVertex> force,
                                          real_t dt) {
  if (force.size() != static_cast<std::size_t>(n_) + 1 ||
      !std::isfinite(dt) || dt < 0) {
    throw std::invalid_argument("invalid force count or coupling timestep");
  }
  std::vector<Vec3> loads(force.size());
  for (std::size_t i = 0; i < force.size(); ++i) {
    const double measure = (i == 0 || i == force.size() - 1) ? h_ / 2.0 : h_;
    loads[i] = measure * Vec3(force[i].x, force[i].y, force[i].z);
    if (!loads[i].allFinite())
      throw std::invalid_argument("force densities must be finite");
  }
  if (dt > 0) {
    state_ = integrateState(state_, dt, loads);
  }
  mesh.GetForces() = force;
  syncMeshes();
}

Eigen::VectorXd BendingRodSegments::integrateState(
    const Eigen::VectorXd &state, double dt, const std::vector<Vec3> &loads) const {
  const double step = integration_step_ > 0 ? integration_step_ : dt;
  return integrateRK45(*this, state, dt, step, relative_tolerance_,
                        absolute_tolerance_, adaptive_, loads);
}

void BendingRodSegments::pack_state(IBModelState &s) const {
  s.ints = {1, n_, static_cast<int64_t>(boundary_), adaptive_ ? 1 : 0};
  s.bytes.clear();
  s.reals = {static_cast<real_t>(length_), static_cast<real_t>(integration_step_),
             static_cast<real_t>(relative_tolerance_), static_cast<real_t>(absolute_tolerance_)};
  for (int i = 0; i < 3; ++i) s.reals.push_back(gravity_[i]);
  for (int i = 0; i < 3; ++i) s.reals.push_back(pin_[i]);
  for (const auto *values : {&masses_, &inertias_, &bending_stiffness_, &bending_damping_})
    s.reals.insert(s.reals.end(), values->begin(), values->end());
  s.reals.insert(s.reals.end(), state_.data(), state_.data() + state_.size());
  for (const auto &force : mesh.GetForces()) {
    s.reals.push_back(force.x);
    s.reals.push_back(force.y);
    s.reals.push_back(force.z);
  }
}

void BendingRodSegments::unpack_state(const IBModelState &s) {
  IBModelState expected;
  BendingRodSegments::pack_state(expected);
  const std::size_t prefix = 10 + masses_.size() + inertias_.size() +
                             bending_stiffness_.size() + bending_damping_.size();
  if (s.ints != expected.ints || s.reals.size() != expected.reals.size() ||
      !s.bytes.empty() ||
      !std::equal(expected.reals.begin(), expected.reals.begin() + prefix, s.reals.begin())) {
    throw std::invalid_argument("incompatible BendingRodSegments checkpoint");
  }
  for (const auto value : s.reals) {
    if (!std::isfinite(value))
      throw std::invalid_argument("non-finite checkpoint state");
  }
  Eigen::VectorXd restored(stateSize());
  std::size_t k = prefix;
  for (int i = 0; i < stateSize(); ++i) restored[i] = s.reals[k++];
  unpack(restored);
  state_ = restored;
  for (auto &force : mesh.GetForces()) {
    force = {s.reals[k], s.reals[k + 1], s.reals[k + 2]};
    k += 3;
  }
  syncMeshes();
}

Eigen::VectorXd BendingRodSegments::initialState(const std::vector<Vec3> &tangents,
                                     const std::vector<Vec3> &angular_velocity,
                                     Vec3 center_of_mass,
                                     Vec3 center_velocity) const {
  const std::vector<Vec3> tangent_values =
      broadcastLinkValues(tangents, Vec3::UnitX(), "tangents");
  const std::vector<Vec3> omega_values =
      broadcastLinkValues(angular_velocity, Vec3::Zero(), "angular_velocity");
  if (!center_of_mass.allFinite() || !center_velocity.allFinite()) {
    throw std::invalid_argument("center position and velocity must be finite");
  }

  Eigen::VectorXd state = Eigen::VectorXd::Zero(stateSize());
  if (isFree()) {
    state.head<3>() = center_of_mass;
  }
  for (int i = 0; i < n_; ++i) {
    if (!tangent_values[i].allFinite() || !omega_values[i].allFinite() ||
        tangent_values[i].norm() == 0.0) {
      throw std::invalid_argument(
          "initial tangents and velocities must be finite");
    }
    const Vec3 tangent = tangent_values[i].normalized();
    state.segment<3>(tangentOffset() + 3 * i) = tangent;
    state.segment<3>(omegaOffset() + 3 * i) =
        omega_values[i] - tangent * tangent.dot(omega_values[i]);
  }
  if (isFree()) {
    state.segment<3>(velocityOffset()) = center_velocity;
  }
  return state;
}

Eigen::VectorXd BendingRodSegments::rhs(double /*time*/,
                            const Eigen::VectorXd &state,
                            const std::vector<Vec3> &nodal_forces) const {
  if (!nodal_forces.empty() && nodal_forces.size() != static_cast<std::size_t>(n_) + 1) {
    throw std::invalid_argument("nodal forces must have links + 1 entries");
  }
  for (const auto &force : nodal_forces) {
    if (!force.allFinite())
      throw std::invalid_argument("nodal forces must be finite");
  }
  const StateParts parts = unpack(state);
  const std::vector<Vec3> alpha =
      angularAcceleration(parts.tangents, parts.omega, nodal_forces);
  Eigen::VectorXd derivative = Eigen::VectorXd::Zero(stateSize());

  if (isFree()) {
    derivative.head<3>() = parts.center_velocity;
    derivative.segment<3>(velocityOffset()) = gravity_;
    for (const auto &force : nodal_forces)
      derivative.segment<3>(velocityOffset()) += force / total_mass_;
  }
  for (int i = 0; i < n_; ++i) {
    derivative.segment<3>(tangentOffset() + 3 * i) =
        parts.omega[i].cross(parts.tangents[i]);
    derivative.segment<3>(omegaOffset() + 3 * i) = alpha[i];
  }
  return derivative;
}

Eigen::MatrixXd BendingRodSegments::nodes(const Eigen::VectorXd &state) const {
  const StateParts parts = unpack(state);
  Eigen::MatrixXd positions(n_ + 1, 3);
  Vec3 start = pin_;
  if (isFree()) {
    start = parts.center;
    for (int i = 0; i < n_; ++i) {
      start -= com_weights_[i] * parts.tangents[i];
    }
  }
  positions.row(0) = start.transpose();
  for (int i = 0; i < n_; ++i) {
    positions.row(i + 1) =
        positions.row(i) + (h_ * parts.tangents[i]).transpose();
  }
  return positions;
}

double BendingRodSegments::energy(const Eigen::VectorXd &state) const {
  const StateParts parts = unpack(state);
  std::vector<Vec3> center_velocity(static_cast<std::size_t>(n_));
  Vec3 endpoint_velocity = Vec3::Zero();
  Vec3 weighted_velocity = Vec3::Zero();
  for (int i = 0; i < n_; ++i) {
    const Vec3 tangent = parts.tangents[i];
    center_velocity[i] =
        endpoint_velocity + parts.omega[i].cross(0.5 * h_ * tangent);
    endpoint_velocity += parts.omega[i].cross(h_ * tangent);
    weighted_velocity += masses_[i] * center_velocity[i];
  }
  if (isFree()) {
    const Vec3 shift = parts.center_velocity - weighted_velocity / total_mass_;
    for (Vec3 &velocity : center_velocity) {
      velocity += shift;
    }
  }

  double kinetic = 0.0;
  for (int i = 0; i < n_; ++i) {
    kinetic += 0.5 * masses_[i] * center_velocity[i].squaredNorm();
    kinetic += 0.5 * inertias_[i] * parts.omega[i].squaredNorm();
  }

  double elastic = 0.0;
  for (int i = 1; i < n_; ++i) {
    const double cosine =
        std::clamp(parts.tangents[i - 1].dot(parts.tangents[i]), -1.0, 1.0);
    const double sine = parts.tangents[i - 1].cross(parts.tangents[i]).norm();
    const double angle = std::atan2(sine, cosine);
    elastic += 0.5 * bending_stiffness_[i - 1] * angle * angle;
  }

  const Eigen::MatrixXd positions = nodes(state);
  double potential = 0.0;
  for (int i = 0; i < n_; ++i) {
    const Vec3 center =
        0.5 * (positions.row(i) + positions.row(i + 1)).transpose();
    potential -= masses_[i] * gravity_.dot(center);
  }
  return kinetic + elastic + potential;
}

void BendingRodSegments::normalizeState(Eigen::VectorXd &state) const {
  if (state.size() != stateSize()) {
    throw std::invalid_argument("state has the wrong size");
  }
  for (int i = 0; i < n_; ++i) {
    Vec3 tangent = state.segment<3>(tangentOffset() + 3 * i);
    const double norm = tangent.norm();
    if (!(norm > 0.0) || !std::isfinite(norm)) {
      throw std::runtime_error("invalid tangent during integration");
    }
    tangent /= norm;
    state.segment<3>(tangentOffset() + 3 * i) = tangent;
    Vec3 omega = state.segment<3>(omegaOffset() + 3 * i);
    state.segment<3>(omegaOffset() + 3 * i) =
        omega - tangent * tangent.dot(omega);
  }
}

void BendingRodSegments::validateNonnegative(double value, const char *name) {
  if (!std::isfinite(value) || value < 0.0) {
    throw std::invalid_argument(std::string(name) +
                                " must be finite and nonnegative");
  }
}

std::vector<Vec3> BendingRodSegments::broadcastLinkValues(const std::vector<Vec3> &supplied,
                                              const Vec3 &fallback,
                                              const char *name) const {
  if (supplied.empty()) {
    return std::vector<Vec3>(static_cast<std::size_t>(n_), fallback);
  }
  if (supplied.size() == 1) {
    return std::vector<Vec3>(static_cast<std::size_t>(n_), supplied.front());
  }
  if (static_cast<int>(supplied.size()) != n_) {
    throw std::invalid_argument(std::string(name) + " must have size 1 or n");
  }
  return supplied;
}

int BendingRodSegments::tangentOffset() const { return isFree() ? 3 : 0; }
int BendingRodSegments::velocityOffset() const { return tangentOffset() + 3 * n_; }
int BendingRodSegments::omegaOffset() const {
  return velocityOffset() + (isFree() ? 3 : 0);
}

BendingRodSegments::StateParts BendingRodSegments::unpack(const Eigen::VectorXd &state) const {
  if (state.size() != stateSize() || !state.allFinite()) {
    throw std::invalid_argument(
        "state has the wrong size or non-finite values");
  }
  StateParts parts;
  parts.tangents.reserve(static_cast<std::size_t>(n_));
  parts.omega.reserve(static_cast<std::size_t>(n_));
  if (isFree()) {
    parts.center = state.head<3>();
    parts.center_velocity = state.segment<3>(velocityOffset());
  }
  for (int i = 0; i < n_; ++i) {
    Vec3 tangent = state.segment<3>(tangentOffset() + 3 * i);
    const double norm = tangent.norm();
    if (!(norm > 0.0)) {
      throw std::invalid_argument("state contains a zero tangent");
    }
    tangent /= norm;
    parts.tangents.push_back(tangent);
    const Vec3 omega = state.segment<3>(omegaOffset() + 3 * i);
    parts.omega.push_back(omega - tangent * tangent.dot(omega));
  }
  return parts;
}

std::vector<Vec3> BendingRodSegments::jointTorques(const std::vector<Vec3> &tangents,
                                       const std::vector<Vec3> &omega) const {
  std::vector<Vec3> torque(static_cast<std::size_t>(n_), Vec3::Zero());
  for (int i = 1; i < n_; ++i) {
    const Vec3 cross = tangents[i - 1].cross(tangents[i]);
    const double sine = cross.norm();
    const double cosine =
        std::clamp(tangents[i - 1].dot(tangents[i]), -1.0, 1.0);
    Vec3 child_torque = Vec3::Zero();
    if (sine > 1e-12) {
      child_torque -=
          bending_stiffness_[i - 1] * std::atan2(sine, cosine) * cross / sine;
    }
    Vec3 mean_tangent = tangents[i - 1] + tangents[i];
    if (mean_tangent.squaredNorm() < 1e-20) {
      mean_tangent = tangents[i - 1];
    }
    mean_tangent.normalize();
    const Mat3 bending_projector =
        Mat3::Identity() - mean_tangent * mean_tangent.transpose();
    child_torque -=
        bending_damping_[i - 1] * bending_projector * (omega[i] - omega[i - 1]);
    torque[i] += child_torque;
    torque[i - 1] -= child_torque;
  }
  return torque;
}

std::vector<Vec3>
BendingRodSegments::angularAcceleration(const std::vector<Vec3> &tangents,
                            const std::vector<Vec3> &omega,
                            const std::vector<Vec3> &nodal_forces) const {
  std::vector<Mat6> inertia(static_cast<std::size_t>(n_), Mat6::Zero());
  std::vector<Vec6> bias(static_cast<std::size_t>(n_), Vec6::Zero());
  std::vector<Mat6> transform(static_cast<std::size_t>(n_), Mat6::Zero());
  std::vector<Vec6> acceleration_bias(static_cast<std::size_t>(n_),
                                      Vec6::Zero());
  std::vector<Mat32> basis(static_cast<std::size_t>(n_));
  std::vector<Mat62> motion(static_cast<std::size_t>(n_), Mat62::Zero());

  for (int i = 0; i < n_; ++i) {
    basis[i] = tangentBasis(tangents[i]);
    motion[i].topRows<3>() = basis[i];
    transform[i].bottomRightCorner<3, 3>() = Mat3::Identity();
    const Vec3 center = 0.5 * h_ * tangents[i];
    const Mat3 center_cross = crossMatrix(center);
    const Mat3 center_inertia =
        inertias_[i] *
        (Mat3::Identity() - tangents[i] * tangents[i].transpose());
    inertia[i].topLeftCorner<3, 3>() =
        center_inertia - masses_[i] * center_cross * center_cross;
    inertia[i].topRightCorner<3, 3>() = masses_[i] * center_cross;
    inertia[i].bottomLeftCorner<3, 3>() = -masses_[i] * center_cross;
    inertia[i].bottomRightCorner<3, 3>() = masses_[i] * Mat3::Identity();

    const Vec3 centripetal = omega[i].cross(omega[i].cross(center));
    const Vec3 effective = masses_[i] * (centripetal - gravity_);
    bias[i].head<3>() =
        omega[i].cross(center_inertia * omega[i]) + center.cross(effective);
    bias[i].tail<3>() = effective;
    if (!nodal_forces.empty()) {
      // Assign endpoint i+1 to the distal end of segment i. Node zero is
      // applied at the proximal end of segment zero. Joint transmission
      // carries these loads to the rest of the articulated chain.
      bias[i].head<3>() -= (h_ * tangents[i]).cross(nodal_forces[i + 1]);
      bias[i].tail<3>() -= nodal_forces[i + 1];
      if (i == 0)
        bias[i].tail<3>() -= nodal_forces[0];
    }
  }
  for (int i = 1; i < n_; ++i) {
    const Vec3 offset = h_ * tangents[i - 1];
    transform[i].bottomLeftCorner<3, 3>() = -crossMatrix(offset);
    acceleration_bias[i].tail<3>() =
        omega[i - 1].cross(omega[i - 1].cross(offset));
  }

  const std::vector<Vec3> torque = jointTorques(tangents, omega);
  std::vector<Mat62> U(static_cast<std::size_t>(n_), Mat62::Zero());
  std::vector<Mat2> d(static_cast<std::size_t>(n_), Mat2::Zero());
  std::vector<Vec2> u(static_cast<std::size_t>(n_), Vec2::Zero());
  const int first_joint = isFree() ? 1 : 0;

  for (int i = n_ - 1; i >= first_joint; --i) {
    U[i] = inertia[i] * motion[i];
    d[i] = motion[i].transpose() * U[i];
    Vec6 net_force = -bias[i];
    net_force.head<3>() += torque[i];
    u[i] = motion[i].transpose() * net_force;
    const Eigen::LDLT<Mat2> factor(d[i]);
    if (factor.info() != Eigen::Success) {
      throw std::runtime_error("singular bending inertia");
    }
    if (i > 0) {
      const Mat6 reduced = inertia[i] - U[i] * factor.solve(U[i].transpose());
      const Vec6 force =
          bias[i] + reduced * acceleration_bias[i] + U[i] * factor.solve(u[i]);
      inertia[i - 1] += transform[i].transpose() * reduced * transform[i];
      bias[i - 1] += transform[i].transpose() * force;
    }
  }

  Vec6 acceleration = Vec6::Zero();
  std::vector<Vec3> alpha(static_cast<std::size_t>(n_), Vec3::Zero());
  if (isFree()) {
    Mat65 base_motion = Mat65::Zero();
    base_motion.topLeftCorner<3, 2>() = basis[0];
    base_motion.bottomRightCorner<3, 3>() = Mat3::Identity();
    const Mat5 base_inertia =
        base_motion.transpose() * inertia[0] * base_motion;
    Vec6 base_force = -bias[0];
    base_force.head<3>() += torque[0];
    const Vec5 base_u = base_motion.transpose() * base_force;
    const Eigen::LDLT<Mat5> factor(base_inertia);
    if (factor.info() != Eigen::Success) {
      throw std::runtime_error("singular floating-base inertia");
    }
    acceleration = base_motion * factor.solve(base_u);
    alpha[0] = acceleration.head<3>();
  }
  for (int i = first_joint; i < n_; ++i) {
    acceleration = transform[i] * acceleration + acceleration_bias[i];
    const Eigen::LDLT<Mat2> factor(d[i]);
    const Vec2 bending_alpha =
        factor.solve(u[i] - U[i].transpose() * acceleration);
    acceleration += motion[i] * bending_alpha;
    alpha[i] = acceleration.head<3>();
  }
  return alpha;
}

Eigen::VectorXd integrateRK45(const BendingRodSegments &beam, Eigen::VectorXd state,
                              double duration, double initial_step,
                              double relative_tolerance,
                              double absolute_tolerance, bool adaptive,
                              const std::vector<Vec3> &nodal_forces) {
  if (!(duration >= 0.0) || !(initial_step > 0.0) ||
      (adaptive && (!(relative_tolerance > 0.0) ||
                    !(absolute_tolerance > 0.0) ||
                    !std::isfinite(relative_tolerance) ||
                    !std::isfinite(absolute_tolerance))) ||
      !std::isfinite(duration) || !std::isfinite(initial_step)) {
    throw std::invalid_argument("invalid RK45 integration parameters");
  }

  if (state.size() != beam.stateSize() || !state.allFinite()) {
    throw std::invalid_argument("state has the wrong size or non-finite values");
  }

  const auto evaluate = [&](double time, const Eigen::VectorXd &value) {
    return beam.rhs(time, value, nodal_forces);
  };

  // Dormand--Prince 5(4). The fifth-order state is accepted and the
  // embedded fourth-order state supplies a local error estimate.
  double time = 0.0;
  double step = std::min(initial_step, duration);
  while (time < duration) {
    step = std::min(step, duration - time);
    if (!(step > 0.0) || time + step == time) {
      throw std::runtime_error("RK45 time step underflow");
    }

    const Eigen::VectorXd k1 = evaluate(time, state);
    const Eigen::VectorXd k2 =
        evaluate(time + step * (1.0 / 5.0), state + step * ((1.0 / 5.0) * k1));
    const Eigen::VectorXd k3 =
        evaluate(time + step * (3.0 / 10.0),
                 state + step * ((3.0 / 40.0) * k1 + (9.0 / 40.0) * k2));
    const Eigen::VectorXd k4 =
        evaluate(time + step * (4.0 / 5.0),
                 state + step * ((44.0 / 45.0) * k1 - (56.0 / 15.0) * k2 +
                                 (32.0 / 9.0) * k3));
    const Eigen::VectorXd k5 = evaluate(
        time + step * (8.0 / 9.0),
        state + step * ((19372.0 / 6561.0) * k1 - (25360.0 / 2187.0) * k2 +
                        (64448.0 / 6561.0) * k3 - (212.0 / 729.0) * k4));
    const Eigen::VectorXd k6 =
        evaluate(time + step,
                 state + step * ((9017.0 / 3168.0) * k1 - (355.0 / 33.0) * k2 +
                                 (46732.0 / 5247.0) * k3 + (49.0 / 176.0) * k4 -
                                 (5103.0 / 18656.0) * k5));
    const Eigen::VectorXd fifth_order =
        state + step * ((35.0 / 384.0) * k1 + (500.0 / 1113.0) * k3 +
                        (125.0 / 192.0) * k4 - (2187.0 / 6784.0) * k5 +
                        (11.0 / 84.0) * k6);
    if (!adaptive) {
      if (!fifth_order.allFinite()) {
        throw std::runtime_error("non-finite state during fixed-step RK45");
      }
      state = fifth_order;
      beam.normalizeState(state);
      time += step;
      continue;
    }

    const Eigen::VectorXd k7 = evaluate(time + step, fifth_order);
    const Eigen::VectorXd fourth_order =
        state + step * ((5179.0 / 57600.0) * k1 + (7571.0 / 16695.0) * k3 +
                        (393.0 / 640.0) * k4 - (92097.0 / 339200.0) * k5 +
                        (187.0 / 2100.0) * k6 + (1.0 / 40.0) * k7);

    const Eigen::ArrayXd scale =
        absolute_tolerance +
        relative_tolerance * state.array().abs().max(fifth_order.array().abs());
    const Eigen::ArrayXd scaled_error =
        (fifth_order - fourth_order).array() / scale;
    const double error = std::sqrt(scaled_error.square().mean());

    if (std::isfinite(error) && error <= 1.0) {
      state = fifth_order;
      beam.normalizeState(state);
      time += step;
    }
    const double factor =
        error == 0.0 ? 5.0 : std::clamp(0.9 * std::pow(error, -0.2), 0.1, 5.0);
    step *= std::isfinite(factor) ? factor : 0.1;
  }
  return state;
}

} // namespace Models
} // namespace ELFF
