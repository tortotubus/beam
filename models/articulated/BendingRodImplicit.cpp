#include "elff/models/articulated/BendingRodImplicit.hpp"

#include <algorithm>
#include <climits>
#include <cmath>
#include <limits>
#include "elff/general/error.hpp"
#include <unsupported/Eigen/AutoDiff>

namespace ELFF {
namespace Models {
namespace {
using Vec2 = Eigen::Vector2d;
using Vec5 = Eigen::Matrix<double, 5, 1>;
using Vec6 = Eigen::Matrix<double, 6, 1>;
using Mat3 = Eigen::Matrix3d;
using Mat5 = Eigen::Matrix<double, 5, 5>;
using Mat6 = Eigen::Matrix<double, 6, 6>;
using Mat62 = Eigen::Matrix<double, 6, 2>;
using Mat65 = Eigen::Matrix<double, 6, 5>;
using Mat2 = Eigen::Matrix2d;
using Mat32 = Eigen::Matrix<double, 3, 2>;

Mat3 crossMatrix(const Vec3 &value) {
  Mat3 result;
  result << 0.0, -value.z(), value.y(), value.z(), 0.0, -value.x(), -value.y(),
      value.x(), 0.0;
  return result;
}

template <typename Value>
std::vector<Value> configuredValues(const std::vector<Value> &supplied,
                                    int count, const Value &fallback,
                                    const char *name) {
  if (supplied.empty()) {
    return std::vector<Value>(static_cast<std::size_t>(count), fallback);
  }
  if (static_cast<int>(supplied.size()) != count) {
    ELFF_ABORT(name << " has the wrong size");
  }
  return supplied;
}

using Basis = Eigen::Matrix<double, 3, 2>;
size_t checkedNodeCount(const BendingRodImplicitConfig &config) {
  if (config.links < 1 || config.links > INT_MAX / 6 - 1 ||
      (config.boundary != Boundary::Pinned &&
       config.boundary != Boundary::Free))
    ELFF_ABORT("invalid implicit rod link count or boundary");
  return static_cast<size_t>(config.links) + 1;
}

Basis tangentBasis(const Vec3 &t) {
  Eigen::Index axis;
  t.cwiseAbs().minCoeff(&axis);
  Basis b;
  b.col(0) = t.cross(Vec3::Unit(axis)).normalized();
  b.col(1) = t.cross(b.col(0));
  return b;
}
template <typename T> using Vector2 = Eigen::Matrix<T, 2, 1>;
template <typename T> using Vector3 = Eigen::Matrix<T, 3, 1>;
template <typename T> struct Chart {
  Vector3<T> tangent, curvature;
  Eigen::Matrix<T, 3, 2> jacobian;
};
template <typename T>
Chart<T> chart(const Vec3 &origin, const Basis &b, const Vector2<T> &u,
               const Vector2<T> &v) {
  using std::cos;
  using std::sin;
  using std::sqrt;
  const T r2 = u.squaredNorm();
  T sinc, cosine, c1, c2;
  // Do not differentiate sqrt(r2) at the chart origin. All these functions
  // are smooth functions of r2 there, including the cosine in the tangent.
  if (r2 < 1e-6) {
    sinc = 1.0 - r2 / 6.0 + r2 * r2 / 120.0;
    cosine = 1.0 - r2 / 2.0 + r2 * r2 / 24.0 - r2 * r2 * r2 / 720.0;
    c1 = -1.0 / 3.0 + r2 / 30.0 - r2 * r2 / 840.0;
    c2 = 1.0 / 15.0 - r2 / 210.0 + r2 * r2 / 7560.0;
  } else {
    const T r = sqrt(r2);
    sinc = sin(r) / r;
    cosine = cos(r);
    c1 = (cosine - sinc) / r2;
    c2 = (-sinc - 3.0 * c1) / r2;
  }
  const Vector3<T> o = origin.template cast<T>();
  const Eigen::Matrix<T, 3, 2> basis = b.template cast<T>();
  const Vector3<T> bu = basis * u;
  const T uv = u.dot(v);
  Chart<T> c;
  c.tangent = cosine * o + sinc * bu;
  c.jacobian = sinc * basis + (-sinc * o + c1 * bu) * u.transpose();
  c.curvature = (-sinc * v.squaredNorm() - c1 * uv * uv) * o +
                2.0 * c1 * uv * (basis * v) +
                (c1 * v.squaredNorm() + c2 * uv * uv) * bu;
  return c;
}

template <typename T>
Vector3<T> jointTorque(const Vector3<T> &a, const Vector3<T> &b,
                       const Vector3<T> &wa, const Vector3<T> &wb,
                       double stiffness, double damping) {
  using std::atan2;
  using std::sqrt;
  const Vector3<T> cross = a.cross(b);
  const T sine2 = cross.squaredNorm();
  const T cosine = a.dot(b);
  Vector3<T> torque = Vector3<T>::Zero();
  if (cosine > 0 && sine2 < 1e-12) {
    // The straight-joint torque has a nonzero derivative, despite its zero
    // value. Use the smooth angle/sin(angle) limit instead of skipping it.
    const T ratio = 1.0 + sine2 / 6.0 + 3.0 * sine2 * sine2 / 40.0;
    torque = -stiffness * ratio * cross;
  } else if (sine2 > 1e-24) {
    const T sine = sqrt(sine2);
    torque = -stiffness * (atan2(sine, cosine) / sine) * cross;
  }
  Vector3<T> mean = a + b;
  if (mean.squaredNorm() < 1e-20)
    mean = a;
  mean.normalize();
  const Vector3<T> relative = wb - wa;
  torque -= damping * (relative - mean * mean.dot(relative));
  return torque;
}
} // namespace

BendingRodImplicit::BendingRodImplicit(const BendingRodImplicitConfig &config)
    : BendingRod(checkedNodeCount(config), config.length), n_(config.links),
      length_(config.length), total_mass_(config.total_mass),
      boundary_(config.boundary), gravity_(config.gravity), pin_(config.pin),
      spectral_radius_(config.spectral_radius),
      tolerance_(config.nonlinear_tolerance),
      max_iterations_(config.max_iterations) {
  if (!std::isfinite(spectral_radius_) || spectral_radius_ < 0 ||
      spectral_radius_ > 1 || !std::isfinite(tolerance_) || tolerance_ <= 0 ||
      max_iterations_ < 1)
    ELFF_ABORT("invalid generalized-alpha solver parameters");
  if (n_ < 1 || !(length_ > 0.0) || !(total_mass_ > 0.0) ||
      !std::isfinite(length_) || !std::isfinite(total_mass_)) {
    ELFF_ABORT(
        "links, length, and total_mass must be positive");
  }
  if (!std::isfinite(config.ei) || config.ei < 0.0 || !gravity_.allFinite() ||
      !pin_.allFinite()) {
    ELFF_ABORT(
        "beam parameters must be finite and nonnegative");
  }
  h_ = length_ / static_cast<double>(n_);

  masses_ = configuredValues(config.masses, n_, total_mass_ / n_, "masses");
  total_mass_ = 0.0;
  for (double mass : masses_) {
    if (!std::isfinite(mass) || mass <= 0.0) {
      ELFF_ABORT("masses must be finite and positive");
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
      ELFF_ABORT(
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

  setState(initialState());
  Eigen::MatrixXd weights(n_, n_);
  for (int body = 0; body < n_; ++body)
    for (int link = 0; link < n_; ++link)
      weights(body, link) = (link < body ? h_ : (link == body ? h_ / 2 : 0.0)) -
                            (isFree() ? com_weights_[link] : 0.0);
  const Eigen::Map<const Eigen::VectorXd> masses(masses_.data(), n_);
  tangent_mass_ = weights.transpose() * masses.asDiagonal() * weights;
  for (int i = 0; i < n_; ++i)
    tangent_mass_(i, i) += inertias_[i];
}

void BendingRodImplicit::setState(const Eigen::VectorXd &state) {
  unpack(state);
  Eigen::VectorXd normalized = state;
  normalizeState(normalized);
  state_.swap(normalized);
  acceleration_.resize(0);
  previous_dt_ = 0;
  update_mesh();
}

Eigen::VectorXd
BendingRodImplicit::integrateState(const Eigen::VectorXd &state, double dt,
                                   const std::vector<Vec3> &loads) const {
  if (!std::isfinite(dt) || dt <= 0)
    ELFF_ABORT(
        "generalized-alpha timestep must be finite and positive");
  const auto parts = unpack(state);
  const auto physical =
      rhs(0, state, loads); // Validate loads, initialize history.
  const int shift = isFree() ? 3 : 0;
  const int size = shift + 2 * n_;
  std::vector<Basis> basis(n_);
  Eigen::VectorXd velocity(size), old_acceleration(size);
  const bool history =
      acceleration_.size() == shift + 3 * n_ && previous_dt_ == dt;
  if (isFree()) {
    velocity.head<3>() = parts.center_velocity;
    old_acceleration.head<3>() =
        history ? Vec3(acceleration_.head<3>())
                : Vec3(physical.segment<3>(velocityOffset()));
  }
  for (int i = 0; i < n_; ++i) {
    basis[i] = tangentBasis(parts.tangents[i]);
    velocity.segment<2>(shift + 2 * i) =
        basis[i].transpose() * parts.omega[i].cross(parts.tangents[i]);
    // At the chart origin its Hessian is normal to the sphere, so its
    // contribution vanishes in these two coordinate components.
    const Vec3 ambient = history
                             ? Vec3(acceleration_.segment<3>(shift + 3 * i))
                             : Vec3(physical.segment<3>(omegaOffset() + 3 * i)
                                        .cross(parts.tangents[i]));
    old_acceleration.segment<2>(shift + 2 * i) = basis[i].transpose() * ambient;
  }

  const double am = (2.0 - spectral_radius_) / (1.0 + spectral_radius_);
  const double af = 1.0 / (1.0 + spectral_radius_);
  const double gamma = 0.5 + am - af;
  const double beta = 0.25 * std::pow(1.0 + am - af, 2);
  const Eigen::VectorXd predictor =
      dt * velocity + dt * dt * (0.5 - beta) * old_acceleration;
  const Eigen::VectorXd velocity_predictor =
      velocity + dt * (1.0 - gamma) * old_acceleration;

  auto makeState = [&](const Eigen::VectorXd &q, const Eigen::VectorXd &v) {
    Eigen::VectorXd result = state;
    if (isFree()) {
      result.head<3>() = parts.center + q.head<3>();
      result.segment<3>(velocityOffset()) = v.head<3>();
    }
    for (int i = 0; i < n_; ++i) {
      const auto c = chart<double>(parts.tangents[i], basis[i],
                                   q.segment<2>(shift + 2 * i),
                                   v.segment<2>(shift + 2 * i));
      result.segment<3>(tangentOffset() + 3 * i) = c.tangent;
      result.segment<3>(omegaOffset() + 3 * i) =
          c.tangent.cross(c.jacobian * v.segment<2>(shift + 2 * i));
    }
    return result;
  };
  // Ambient tangent coordinates have a constant scalar mass matrix. Form
  // the momentum residual directly, avoiding an inverse mass inside Newton.
  std::vector<Vec3> external(n_, Vec3::Zero());
  for (int i = 0; i < n_; ++i) {
    for (int body = 0; body < n_; ++body)
      external[i] += masses_[body] * gravity_ *
                     ((i < body ? h_ : (i == body ? h_ / 2 : 0.0)) -
                      (isFree() ? com_weights_[i] : 0.0));
    for (int node = 0; node < static_cast<int>(loads.size()); ++node)
      external[i] +=
          ((node > i ? h_ : 0.0) - (isFree() ? com_weights_[i] : 0.0)) *
          loads[node];
  }
  // Unknown x=dt^2*a_new keeps perturbations in position units.
  auto residual = [&](const Eigen::VectorXd &x) -> Eigen::VectorXd {
    const Eigen::VectorXd q = af * (predictor + beta * x);
    const Eigen::VectorXd v =
        (1.0 - af) * velocity + af * (velocity_predictor + gamma / dt * x);
    const Eigen::VectorXd weighted =
        am * x + (1.0 - am) * dt * dt * old_acceleration;
    std::vector<Chart<double>> charts(n_);
    std::vector<Vec3> tangents(n_), omega(n_);
    Eigen::MatrixXd ambient(n_, 3);
    for (int i = 0; i < n_; ++i) {
      charts[i] = chart<double>(parts.tangents[i], basis[i],
                                q.segment<2>(shift + 2 * i),
                                v.segment<2>(shift + 2 * i));
      const auto &c = charts[i];
      tangents[i] = c.tangent;
      omega[i] = c.tangent.cross(c.jacobian * v.segment<2>(shift + 2 * i));
      ambient.row(i) = (c.jacobian * weighted.segment<2>(shift + 2 * i) +
                        dt * dt * c.curvature)
                           .transpose();
    }
    std::vector<Vec3> torque(n_, Vec3::Zero());
    for (int i = 1; i < n_; ++i) {
      const Vec3 joint = jointTorque<double>(
          tangents[i - 1], tangents[i], omega[i - 1], omega[i],
          bending_stiffness_[i - 1], bending_damping_[i - 1]);
      torque[i] += joint;
      torque[i - 1] -= joint;
    }
    const Eigen::MatrixXd momentum = tangent_mass_ * ambient;
    Eigen::VectorXd result(size);
    if (isFree())
      result.head<3>() =
          weighted.head<3>() - dt * dt * physical.segment<3>(velocityOffset());
    for (int i = 0; i < n_; ++i)
      result.segment<2>(shift + 2 * i) =
          charts[i].jacobian.transpose() *
          (Vec3(momentum.row(i).transpose()) -
           dt * dt * (torque[i].cross(tangents[i]) + external[i])) /
          tangent_mass_(i, i);
    return result;
  };

  // Start from constant-velocity positions

  Eigen::VectorXd x = -(0.5 - beta) / beta * dt * dt * old_acceleration;
  Eigen::VectorXd r = residual(x);
  const double threshold =
      tolerance_ * dt * (1.0 + velocity.lpNorm<Eigen::Infinity>());
  bool converged = false;
  for (int iteration = 0; iteration <= max_iterations_; ++iteration) {
    if (!r.allFinite())
      break;
    if (r.lpNorm<Eigen::Infinity>() <= threshold) {
      converged = true;
      break;
    }
    if (iteration == max_iterations_)
      break;
    using AD = Eigen::AutoDiffScalar<Eigen::Matrix<double, 4, 1>>;
    std::vector<Chart<AD>> left(n_), right(n_);
    std::vector<Vector3<AD>> ambient(n_), omega_left(n_), omega_right(n_);
    for (int i = 0; i < n_; ++i) {
      for (int side = 0; side < 2; ++side) {
        Vector2<AD> q, v, weighted;
        for (int k = 0; k < 2; ++k) {
          const int index = shift + 2 * i + k;
          AD variable(x[index]);
          variable.derivatives().setZero();
          variable.derivatives()[2 * side + k] = 1.0;
          q[k] = af * (predictor[index] + beta * variable);
          v[k] = (1.0 - af) * velocity[index] +
                 af * (velocity_predictor[index] + gamma / dt * variable);
          weighted[k] =
              am * variable + (1.0 - am) * dt * dt * old_acceleration[index];
        }
        auto &c = side == 0 ? left[i] : right[i];
        c = chart<AD>(parts.tangents[i], basis[i], q, v);
        (side == 0 ? omega_left[i] : omega_right[i]) =
            c.tangent.cross(c.jacobian * v);
        if (side == 1)
          ambient[i] = c.jacobian * weighted + dt * dt * c.curvature;
      }
    }
    Eigen::MatrixXd jacobian = Eigen::MatrixXd::Zero(size, size);
    if (isFree())
      jacobian.topLeftCorner<3, 3>().diagonal().setConstant(am);
    auto assemble = [&](int row, int first, int second,
                        const Vector2<AD> &value) {
      for (int k = 0; k < 2; ++k)
        for (int d = 0; d < 2; ++d) {
          jacobian(shift + 2 * row + k, shift + 2 * first + d) +=
              value[k].derivatives()[d] / tangent_mass_(row, row);
          jacobian(shift + 2 * row + k, shift + 2 * second + d) +=
              value[k].derivatives()[2 + d] / tangent_mass_(row, row);
        }
    };
    for (int i = 0; i < n_; ++i) {
      for (int j = 0; j < n_; ++j) {
        const Vector2<AD> inertia =
            left[i].jacobian.transpose() * (tangent_mass_(i, j) * ambient[j]);
        assemble(i, i, j, inertia);
      }
      const Vector2<AD> force =
          -dt * dt * left[i].jacobian.transpose() * external[i].cast<AD>();
      assemble(i, i, i, force);
    }
    for (int i = 1; i < n_; ++i) {
      const auto &a = left[i - 1];
      const auto &b = right[i];
      const Vector3<AD> torque = jointTorque<AD>(
          a.tangent, b.tangent, omega_left[i - 1], omega_right[i],
          bending_stiffness_[i - 1], bending_damping_[i - 1]);
      const Vector2<AD> parent =
          dt * dt * a.jacobian.transpose() * torque.cross(a.tangent);
      const Vector2<AD> child =
          -dt * dt * b.jacobian.transpose() * torque.cross(b.tangent);
      assemble(i - 1, i - 1, i, parent);
      assemble(i, i - 1, i, child);
    }
    const Eigen::VectorXd delta = jacobian.partialPivLu().solve(-r);
    if (!delta.allFinite())
      break;
    bool accepted = false;
    for (double scale = 1.0; scale >= 1.0 / 1024; scale *= 0.5) {
      const Eigen::VectorXd candidate = x + scale * delta;
      const Eigen::VectorXd next = residual(candidate);
      if (next.allFinite() && next.norm() < r.norm()) {
        x = candidate;
        r = next;
        accepted = true;
        break;
      }
    }
    if (!accepted)
      break;
  }
  if (!converged) {
    ELFF_ABORT("generalized-alpha Newton solve failed (residual "
            << r.lpNorm<Eigen::Infinity>() << ", tolerance " << threshold
            << "); reduce timestep");
  }

  const Eigen::VectorXd q = predictor + beta * x;
  const Eigen::VectorXd v = velocity_predictor + gamma / dt * x;
  Eigen::VectorXd result = makeState(q, v);
  Eigen::VectorXd next_acceleration(shift + 3 * n_);
  if (isFree())
    next_acceleration.head<3>() = x.head<3>() / (dt * dt);
  for (int i = 0; i < n_; ++i) {
    const auto c =
        chart<double>(parts.tangents[i], basis[i], q.segment<2>(shift + 2 * i),
                      v.segment<2>(shift + 2 * i));
    next_acceleration.segment<3>(shift + 3 * i) =
        c.jacobian * (x.segment<2>(shift + 2 * i) / (dt * dt)) + c.curvature;
  }
  if (!result.allFinite() || !next_acceleration.allFinite())
    ELFF_ABORT("non-finite generalized-alpha state");
  acceleration_ = next_acceleration;
  previous_dt_ = dt;
  return result;
}

void BendingRodImplicit::advance(double dt, const std::vector<Vec3> &loads) {
  if (!std::isfinite(dt) || dt < 0 ||
      (!loads.empty() && loads.size() != static_cast<size_t>(n_) + 1))
    ELFF_ABORT("invalid timestep or endpoint force count");
  for (const auto &force : loads)
    if (!force.allFinite())
      ELFF_ABORT("endpoint forces must be finite");
  if (dt == 0)
    return;
  Eigen::VectorXd next = integrateState(state_, dt, loads);
  state_.swap(next);
  update_mesh();
}

void BendingRodImplicit::pack_state(IBModelState &s) const {
  s.ints = {1, n_, static_cast<int64_t>(boundary_), max_iterations_,
            acceleration_.size() != 0 ? 1 : 0};
  s.bytes.clear();
  s.reals = {length_, spectral_radius_, tolerance_};
  for (int i = 0; i < 3; ++i)
    s.reals.push_back(gravity_[i]);
  for (int i = 0; i < 3; ++i)
    s.reals.push_back(pin_[i]);
  for (const auto *values :
       {&masses_, &inertias_, &bending_stiffness_, &bending_damping_})
    s.reals.insert(s.reals.end(), values->begin(), values->end());
  s.reals.insert(s.reals.end(), state_.data(), state_.data() + state_.size());
  s.reals.push_back(previous_dt_);
  if (acceleration_.size())
    s.reals.insert(s.reals.end(), acceleration_.data(),
                   acceleration_.data() + acceleration_.size());
}

void BendingRodImplicit::unpack_state(const IBModelState &s) {
  IBModelState expected;
  BendingRodImplicit::pack_state(expected);
  const size_t prefix = 9 + masses_.size() + inertias_.size() +
                        bending_stiffness_.size() + bending_damping_.size();
  if (s.ints.size() != 5 ||
      !std::equal(s.ints.begin(), s.ints.begin() + 4, expected.ints.begin()) ||
      (s.ints[4] != 0 && s.ints[4] != 1) || !s.bytes.empty())
    ELFF_ABORT("incompatible generalized-alpha checkpoint");
  const int count = s.ints[4] ? 3 * n_ + (isFree() ? 3 : 0) : 0;
  if (s.reals.size() != prefix + stateSize() + 1 + count ||
      !std::equal(expected.reals.begin(), expected.reals.begin() + prefix,
                  s.reals.begin()))
    ELFF_ABORT(
        "incompatible generalized-alpha checkpoint parameters");
  for (auto value : s.reals)
    if (!std::isfinite(value))
      ELFF_ABORT("non-finite generalized-alpha checkpoint");
  const double dt = s.reals[prefix + stateSize()];
  if ((count && dt <= 0) || (!count && dt != 0))
    ELFF_ABORT("invalid generalized-alpha checkpoint history");
  Eigen::VectorXd restored(stateSize()), acceleration(count);
  for (int i = 0; i < stateSize(); ++i)
    restored[i] = s.reals[prefix + i];
  unpack(restored);
  for (int i = 0; i < count; ++i)
    acceleration[i] = s.reals[prefix + stateSize() + 1 + i];
  state_.swap(restored);
  acceleration_.swap(acceleration);
  previous_dt_ = dt;
  update_mesh();
}

int BendingRodImplicit::stateSize() const {
  return 6 * n_ + (isFree() ? 6 : 0);
}
bool BendingRodImplicit::isFree() const { return boundary_ == Boundary::Free; }

Eigen::VectorXd
BendingRodImplicit::initialState(const std::vector<Vec3> &tangents,
                                 const std::vector<Vec3> &angular_velocity,
                                 Vec3 center_of_mass,
                                 Vec3 center_velocity) const {
  const std::vector<Vec3> tangent_values =
      broadcastLinkValues(tangents, Vec3::UnitX(), "tangents");
  const std::vector<Vec3> omega_values =
      broadcastLinkValues(angular_velocity, Vec3::Zero(), "angular_velocity");
  if (!center_of_mass.allFinite() || !center_velocity.allFinite()) {
    ELFF_ABORT("center position and velocity must be finite");
  }

  Eigen::VectorXd state = Eigen::VectorXd::Zero(stateSize());
  if (isFree()) {
    state.head<3>() = center_of_mass;
  }
  for (int i = 0; i < n_; ++i) {
    if (!tangent_values[i].allFinite() || !omega_values[i].allFinite() ||
        tangent_values[i].norm() == 0.0) {
      ELFF_ABORT(
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

Eigen::VectorXd
BendingRodImplicit::rhs(double /*time*/, const Eigen::VectorXd &state,
                        const std::vector<Vec3> &nodal_forces) const {
  if (!nodal_forces.empty() &&
      nodal_forces.size() != static_cast<std::size_t>(n_) + 1) {
    ELFF_ABORT("nodal forces must have links + 1 entries");
  }
  for (const auto &force : nodal_forces) {
    if (!force.allFinite())
      ELFF_ABORT("nodal forces must be finite");
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

Eigen::MatrixXd BendingRodImplicit::nodes(const Eigen::VectorXd &state) const {
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

double BendingRodImplicit::energy(const Eigen::VectorXd &state) const {
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

void BendingRodImplicit::normalizeState(Eigen::VectorXd &state) const {
  if (state.size() != stateSize()) {
    ELFF_ABORT("state has the wrong size");
  }
  for (int i = 0; i < n_; ++i) {
    Vec3 tangent = state.segment<3>(tangentOffset() + 3 * i);
    const double norm = tangent.norm();
    if (!(norm > 0.0) || !std::isfinite(norm)) {
      ELFF_ABORT("invalid tangent during integration");
    }
    tangent /= norm;
    state.segment<3>(tangentOffset() + 3 * i) = tangent;
    Vec3 omega = state.segment<3>(omegaOffset() + 3 * i);
    state.segment<3>(omegaOffset() + 3 * i) =
        omega - tangent * tangent.dot(omega);
  }
}

void BendingRodImplicit::validateNonnegative(double value, const char *name) {
  if (!std::isfinite(value) || value < 0.0) {
    ELFF_ABORT(name << " must be finite and nonnegative");
  }
}

std::vector<Vec3>
BendingRodImplicit::broadcastLinkValues(const std::vector<Vec3> &supplied,
                                        const Vec3 &fallback,
                                        const char *name) const {
  if (supplied.empty()) {
    return std::vector<Vec3>(static_cast<std::size_t>(n_), fallback);
  }
  if (supplied.size() == 1) {
    return std::vector<Vec3>(static_cast<std::size_t>(n_), supplied.front());
  }
  if (static_cast<int>(supplied.size()) != n_) {
    ELFF_ABORT(name << " must have size 1 or n");
  }
  return supplied;
}

int BendingRodImplicit::tangentOffset() const { return isFree() ? 3 : 0; }
int BendingRodImplicit::velocityOffset() const {
  return tangentOffset() + 3 * n_;
}
int BendingRodImplicit::omegaOffset() const {
  return velocityOffset() + (isFree() ? 3 : 0);
}

BendingRodImplicit::StateParts
BendingRodImplicit::unpack(const Eigen::VectorXd &state) const {
  if (state.size() != stateSize() || !state.allFinite()) {
    ELFF_ABORT(
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
      ELFF_ABORT("state contains a zero tangent");
    }
    tangent /= norm;
    parts.tangents.push_back(tangent);
    const Vec3 omega = state.segment<3>(omegaOffset() + 3 * i);
    parts.omega.push_back(omega - tangent * tangent.dot(omega));
  }
  return parts;
}

std::vector<Vec3>
BendingRodImplicit::jointTorques(const std::vector<Vec3> &tangents,
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

std::vector<Vec3> BendingRodImplicit::angularAcceleration(
    const std::vector<Vec3> &tangents, const std::vector<Vec3> &omega,
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
      ELFF_ABORT("singular bending inertia");
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
      ELFF_ABORT("singular floating-base inertia");
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

Eigen::MatrixXd
BendingRodImplicit::velocities(const Eigen::VectorXd &state) const {
  const auto parts = unpack(state);
  std::vector<Vec3> velocities(static_cast<std::size_t>(n_) + 1, Vec3::Zero());
  Vec3 weighted_velocity = Vec3::Zero();
  for (int i = 0; i < n_; ++i) {
    velocities[i + 1] =
        velocities[i] + parts.omega[i].cross(h_ * parts.tangents[i]);
    weighted_velocity += 0.5 * masses_[i] * (velocities[i] + velocities[i + 1]);
  }
  const Vec3 shift =
      isFree() ? Vec3(parts.center_velocity - weighted_velocity / total_mass_)
               : Vec3::Zero();
  Eigen::MatrixXd result(n_ + 1, 3);
  for (int i = 0; i <= n_; ++i)
    result.row(i) = (velocities[i] + shift).transpose();
  return result;
}

void BendingRodImplicit::update_mesh() {
  const auto positions = nodes(state_);
  const auto velocity = velocities(state_);
  for (size_t i = 0; i < mesh_.get_nodes(); ++i)
    for (size_t d = 0; d < 3; ++d) {
      mesh_.get_centerline(i)[d] = static_cast<real_t>(positions(i, d));
      mesh_.get_centerline_velocity(i)[d] = static_cast<real_t>(velocity(i, d));
    }
}

void BendingRodImplicit::solve(real_t dt, const std::vector<Vector> &density) {
  if (density.size() != static_cast<size_t>(n_) + 1)
    ELFF_ABORT(
        "force density count must match rod mesh nodes");
  std::vector<Vec3> forces(density.size());
  for (size_t i = 0; i < density.size(); ++i) {
    const double weight = (i == 0 || i + 1 == density.size()) ? h_ / 2 : h_;
    forces[i] = weight * Vec3(density[i][0], density[i][1], density[i][2]);
  }
  advance(dt, forces);
}

void BendingRodImplicit::apply_initial_condition(
    const BendingRodMesh &initial) {
  if (initial.get_nodes() != static_cast<size_t>(n_) + 1 ||
      initial.get_length() != mesh_.get_length() ||
      initial.get_centerline().size() != initial.get_nodes() ||
      initial.get_centerline_velocity().size() != initial.get_nodes())
    ELFF_ABORT(
        "initial rod mesh discretization does not match solver");
  std::vector<Vec3> positions(n_ + 1), velocity(n_ + 1), tangents(n_),
      omega(n_);
  double position_scale = length_, velocity_scale = 1;
  for (int i = 0; i <= n_; ++i) {
    const auto &x = initial.get_centerline(i);
    const auto &v = initial.get_centerline_velocity(i);
    positions[i] = Vec3(x[0], x[1], x[2]);
    velocity[i] = Vec3(v[0], v[1], v[2]);
    if (!positions[i].allFinite() || !velocity[i].allFinite())
      ELFF_ABORT("initial rod mesh must be finite");
    position_scale =
        std::max(position_scale, positions[i].cwiseAbs().maxCoeff());
    velocity_scale =
        std::max(velocity_scale, velocity[i].cwiseAbs().maxCoeff());
  }
  const double tolerance = 64 * std::numeric_limits<real_t>::epsilon();
  if (!isFree() && ((positions[0] - pin_).norm() > tolerance * position_scale ||
                    velocity[0].norm() > tolerance * velocity_scale))
    ELFF_ABORT("initial rod mesh must respect the fixed pin");
  Vec3 center = Vec3::Zero(), center_velocity = Vec3::Zero();
  for (int i = 0; i < n_; ++i) {
    const Vec3 edge = positions[i + 1] - positions[i];
    const double length = edge.norm();
    if (!std::isfinite(length) || length <= 0 ||
        std::abs(length - h_) > tolerance * position_scale)
      ELFF_ABORT(
          "initial rod mesh must have inextensible equal-length segments");
    tangents[i] = edge / length;
    const Vec3 delta_velocity = velocity[i + 1] - velocity[i];
    if (std::abs(tangents[i].dot(delta_velocity)) > tolerance * velocity_scale)
      ELFF_ABORT(
          "initial rod velocities must preserve segment lengths");
    omega[i] = tangents[i].cross(delta_velocity / h_);
    center += masses_[i] * (0.5 * positions[i] + 0.5 * positions[i + 1]);
    center_velocity += masses_[i] * (0.5 * velocity[i] + 0.5 * velocity[i + 1]);
  }
  setState(initialState(tangents, omega, center / total_mass_,
                        center_velocity / total_mass_));
}

} // namespace Models
} // namespace ELFF
