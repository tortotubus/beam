#include "elff/models/articulated/BendingRodSegmentsImplicit.hpp"

#include <cmath>
#include <unsupported/Eigen/AutoDiff>
#include <stdexcept>
#include <sstream>

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
Basis tangentBasis(const Vec3 &t)
{
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
Chart<T> chart(const Vec3 &origin, const Basis &b, const Vector2<T> &u, const Vector2<T> &v)
{
  using std::sqrt;
  using std::sin;
  using std::cos;
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
  c.curvature = (-sinc * v.squaredNorm() - c1 * uv * uv) * o
                + 2.0 * c1 * uv * (basis * v)
                + (c1 * v.squaredNorm() + c2 * uv * uv) * bu;
  return c;
}

template <typename T>
Vector3<T> jointTorque(const Vector3<T> &a, const Vector3<T> &b,
                      const Vector3<T> &wa, const Vector3<T> &wb,
                      double stiffness, double damping)
{
  using std::sqrt;
  using std::atan2;
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
  if (mean.squaredNorm() < 1e-20) mean = a;
  mean.normalize();
  const Vector3<T> relative = wb - wa;
  torque -= damping * (relative - mean * mean.dot(relative));
  return torque;
}
} // namespace

BendingRodSegmentsImplicit::BendingRodSegmentsImplicit(
    const BendingRodSegmentsImplicitConfig &config)
  : BendingRodSegments(sharedConfig(config)), spectral_radius_(config.spectral_radius),
    tolerance_(config.nonlinear_tolerance), max_iterations_(config.max_iterations)
{
  if (!std::isfinite(spectral_radius_) || spectral_radius_ < 0 || spectral_radius_ > 1 ||
      !std::isfinite(tolerance_) || tolerance_ <= 0 || max_iterations_ < 1)
    throw std::invalid_argument("invalid generalized-alpha solver parameters");
  Eigen::MatrixXd weights(n_, n_);
  for (int body = 0; body < n_; ++body)
    for (int link = 0; link < n_; ++link)
      weights(body, link) = (link < body ? h_ : (link == body ? h_ / 2 : 0.0))
                           - (isFree() ? com_weights_[link] : 0.0);
  const Eigen::Map<const Eigen::VectorXd> masses(masses_.data(), n_);
  tangent_mass_ = weights.transpose() * masses.asDiagonal() * weights;
  for (int i = 0; i < n_; ++i) tangent_mass_(i, i) += inertias_[i];
}

void BendingRodSegmentsImplicit::setState(const Eigen::VectorXd &state)
{
  BendingRodSegments::setState(state);
  acceleration_.resize(0);
  previous_dt_ = 0;
}

Eigen::VectorXd BendingRodSegmentsImplicit::integrateState(
    const Eigen::VectorXd &state, double dt, const std::vector<Vec3> &loads) const
{
  if (!std::isfinite(dt) || dt <= 0)
    throw std::invalid_argument("generalized-alpha timestep must be finite and positive");
  const auto parts = unpack(state);
  const auto physical = rhs(0, state, loads); // Validate loads, initialize history.
  const int shift = isFree() ? 3 : 0;
  const int size = shift + 2 * n_;
  std::vector<Basis> basis(n_);
  Eigen::VectorXd velocity(size), old_acceleration(size);
  const bool history = acceleration_.size() == shift + 3 * n_ && previous_dt_ == dt;
  if (isFree()) {
    velocity.head<3>() = parts.center_velocity;
    old_acceleration.head<3>() = history ? Vec3(acceleration_.head<3>())
                                        : Vec3(physical.segment<3>(velocityOffset()));
  }
  for (int i = 0; i < n_; ++i) {
    basis[i] = tangentBasis(parts.tangents[i]);
    velocity.segment<2>(shift + 2 * i) = basis[i].transpose() *
        parts.omega[i].cross(parts.tangents[i]);
    // At the chart origin its Hessian is normal to the sphere, so its
    // contribution vanishes in these two coordinate components.
    const Vec3 ambient = history ? Vec3(acceleration_.segment<3>(shift + 3 * i))
        : Vec3(physical.segment<3>(omegaOffset() + 3 * i).cross(parts.tangents[i]));
    old_acceleration.segment<2>(shift + 2 * i) = basis[i].transpose() * ambient;
  }

  // Weights multiply NEW values (PETSc TSALPHA2 convention).
  // https://petsc.org/release/manualpages/TS/TSAlpha2SetRadius/
  const double am = (2.0 - spectral_radius_) / (1.0 + spectral_radius_);
  const double af = 1.0 / (1.0 + spectral_radius_);
  const double gamma = 0.5 + am - af;
  const double beta = 0.25 * std::pow(1.0 + am - af, 2);
  const Eigen::VectorXd predictor = dt * velocity + dt * dt * (0.5 - beta) * old_acceleration;
  const Eigen::VectorXd velocity_predictor = velocity + dt * (1.0 - gamma) * old_acceleration;

  auto makeState = [&](const Eigen::VectorXd &q, const Eigen::VectorXd &v) {
    Eigen::VectorXd result = state;
    if (isFree()) {
      result.head<3>() = parts.center + q.head<3>();
      result.segment<3>(velocityOffset()) = v.head<3>();
    }
    for (int i = 0; i < n_; ++i) {
      const auto c = chart<double>(parts.tangents[i], basis[i], q.segment<2>(shift + 2 * i),
                           v.segment<2>(shift + 2 * i));
      result.segment<3>(tangentOffset() + 3 * i) = c.tangent;
      result.segment<3>(omegaOffset() + 3 * i) = c.tangent.cross(
          c.jacobian * v.segment<2>(shift + 2 * i));
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
      external[i] += ((node > i ? h_ : 0.0) -
                      (isFree() ? com_weights_[i] : 0.0)) * loads[node];
  }
  // Unknown x=dt^2*a_new keeps perturbations in position units.
  auto residual = [&](const Eigen::VectorXd &x) -> Eigen::VectorXd {
    const Eigen::VectorXd q = af * (predictor + beta * x);
    const Eigen::VectorXd v = (1.0 - af) * velocity + af * (velocity_predictor + gamma / dt * x);
    const Eigen::VectorXd weighted = am * x + (1.0 - am) * dt * dt * old_acceleration;
    std::vector<Chart<double>> charts(n_);
    std::vector<Vec3> tangents(n_), omega(n_);
    Eigen::MatrixXd ambient(n_, 3);
    for (int i = 0; i < n_; ++i) {
      charts[i] = chart<double>(parts.tangents[i], basis[i], q.segment<2>(shift + 2 * i),
                         v.segment<2>(shift + 2 * i));
      const auto &c = charts[i];
      tangents[i] = c.tangent;
      omega[i] = c.tangent.cross(c.jacobian * v.segment<2>(shift + 2 * i));
      ambient.row(i) = (c.jacobian * weighted.segment<2>(shift + 2 * i) +
                        dt * dt * c.curvature).transpose();
    }
    std::vector<Vec3> torque(n_, Vec3::Zero());
    for (int i = 1; i < n_; ++i) {
      const Vec3 joint = jointTorque<double>(tangents[i - 1], tangents[i],
          omega[i - 1], omega[i], bending_stiffness_[i - 1], bending_damping_[i - 1]);
      torque[i] += joint;
      torque[i - 1] -= joint;
    }
    const Eigen::MatrixXd momentum = tangent_mass_ * ambient;
    Eigen::VectorXd result(size);
    if (isFree()) result.head<3>() = weighted.head<3>() -
        dt * dt * physical.segment<3>(velocityOffset());
    for (int i = 0; i < n_; ++i)
      result.segment<2>(shift + 2 * i) = charts[i].jacobian.transpose() *
          (Vec3(momentum.row(i).transpose()) - dt * dt *
           (torque[i].cross(tangents[i]) + external[i])) / tangent_mass_(i, i);
    return result;
  };

  // Start from constant-velocity positions, rather than an explicit elastic
  // acceleration predictor that can be enormous for unresolved stiff modes.
  Eigen::VectorXd x = -(0.5 - beta) / beta * dt * dt * old_acceleration;
  Eigen::VectorXd r = residual(x);
  const double threshold = tolerance_ * dt * (1.0 + velocity.lpNorm<Eigen::Infinity>());
  bool converged = false;
  for (int iteration = 0; iteration <= max_iterations_; ++iteration) {
    if (!r.allFinite()) break;
    if (r.lpNorm<Eigen::Infinity>() <= threshold) { converged = true; break; }
    if (iteration == max_iterations_) break;
    // As in the beam element assemblers, seed only the local DoFs. Each
    // pair contribution depends on two two-coordinate segment blocks.
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
          weighted[k] = am * variable + (1.0 - am) * dt * dt * old_acceleration[index];
        }
        auto &c = side == 0 ? left[i] : right[i];
        c = chart<AD>(parts.tangents[i], basis[i], q, v);
        (side == 0 ? omega_left[i] : omega_right[i]) =
            c.tangent.cross(c.jacobian * v);
        if (side == 1) ambient[i] = c.jacobian * weighted + dt * dt * c.curvature;
      }
    }
    Eigen::MatrixXd jacobian = Eigen::MatrixXd::Zero(size, size);
    if (isFree()) jacobian.topLeftCorner<3, 3>().diagonal().setConstant(am);
    auto assemble = [&](int row, int first, int second, const Vector2<AD> &value) {
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
        const Vector2<AD> inertia = left[i].jacobian.transpose() *
                                    (tangent_mass_(i, j) * ambient[j]);
        assemble(i, i, j, inertia);
      }
      const Vector2<AD> force = -dt * dt * left[i].jacobian.transpose() *
                                external[i].cast<AD>();
      assemble(i, i, i, force);
    }
    for (int i = 1; i < n_; ++i) {
      const auto &a = left[i - 1];
      const auto &b = right[i];
      const Vector3<AD> torque = jointTorque<AD>(a.tangent, b.tangent,
          omega_left[i - 1], omega_right[i], bending_stiffness_[i - 1], bending_damping_[i - 1]);
      const Vector2<AD> parent = dt * dt * a.jacobian.transpose() * torque.cross(a.tangent);
      const Vector2<AD> child = -dt * dt * b.jacobian.transpose() * torque.cross(b.tangent);
      assemble(i - 1, i - 1, i, parent);
      assemble(i, i - 1, i, child);
    }
    const Eigen::VectorXd delta = jacobian.partialPivLu().solve(-r);
    if (!delta.allFinite()) break;
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
    if (!accepted) break;
  }
  if (!converged) {
    std::ostringstream message;
    message << "generalized-alpha Newton solve failed (residual " << r.lpNorm<Eigen::Infinity>()
            << ", tolerance " << threshold << "); reduce timestep";
    throw std::runtime_error(message.str());
  }

  const Eigen::VectorXd q = predictor + beta * x;
  const Eigen::VectorXd v = velocity_predictor + gamma / dt * x;
  Eigen::VectorXd result = makeState(q, v);
  Eigen::VectorXd next_acceleration(shift + 3 * n_);
  if (isFree()) next_acceleration.head<3>() = x.head<3>() / (dt * dt);
  for (int i = 0; i < n_; ++i) {
    const auto c = chart<double>(parts.tangents[i], basis[i], q.segment<2>(shift + 2 * i),
                         v.segment<2>(shift + 2 * i));
    next_acceleration.segment<3>(shift + 3 * i) =
        c.jacobian * (x.segment<2>(shift + 2 * i) / (dt * dt)) + c.curvature;
  }
  if (!result.allFinite() || !next_acceleration.allFinite())
    throw std::runtime_error("non-finite generalized-alpha state");
  acceleration_ = next_acceleration;
  previous_dt_ = dt;
  return result;
}

void BendingRodSegmentsImplicit::pack_state(IBModelState &s) const
{
  BendingRodSegments::pack_state(s);
  s.ints.insert(s.ints.end(), {3, max_iterations_, acceleration_.size() != 0 ? 1 : 0});
  s.reals.insert(s.reals.end(), {spectral_radius_, tolerance_, previous_dt_});
  for (int i = 0; i < acceleration_.size(); ++i) s.reals.push_back(acceleration_[i]);
}

void BendingRodSegmentsImplicit::unpack_state(const IBModelState &s)
{
  IBModelState common;
  BendingRodSegments::pack_state(common);
  const auto base_size = common.reals.size();
  if (s.ints.size() != 7 || s.ints[4] != 3 || s.ints[5] != max_iterations_ ||
      (s.ints[6] != 0 && s.ints[6] != 1))
    throw std::invalid_argument("incompatible generalized-alpha checkpoint");
  const int count = s.ints[6] ? 3 * n_ + (isFree() ? 3 : 0) : 0;
  if (s.reals.size() != base_size + 3 + count ||
      s.reals[base_size] != spectral_radius_ || s.reals[base_size + 1] != tolerance_)
    throw std::invalid_argument("incompatible generalized-alpha checkpoint parameters");
  for (const auto value : s.reals)
    if (!std::isfinite(value)) throw std::invalid_argument("non-finite generalized-alpha checkpoint");
  const double dt = s.reals[base_size + 2];
  if ((count && dt <= 0) || (!count && dt != 0))
    throw std::invalid_argument("invalid generalized-alpha checkpoint history");
  Eigen::VectorXd restored(count);
  for (int i = 0; i < count; ++i) restored[i] = s.reals[base_size + 3 + i];
  common = s;
  common.ints.resize(4);
  common.reals.resize(base_size);
  BendingRodSegments::unpack_state(common);
  acceleration_ = restored;
  previous_dt_ = dt;
}

} // namespace Models
} // namespace ELFF
