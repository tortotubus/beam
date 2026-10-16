#pragma once

#include "elff/models/ibm/IBForceCoupled.hpp"

#include <Eigen/Dense>
#include <vector>

namespace ELFF {
namespace Models {

using Vec2 = Eigen::Vector2d;
using Vec3 = Eigen::Vector3d;
using Vec5 = Eigen::Matrix<double, 5, 1>;
using Vec6 = Eigen::Matrix<double, 6, 1>;
using Mat3 = Eigen::Matrix3d;
using Mat5 = Eigen::Matrix<double, 5, 5>;
using Mat6 = Eigen::Matrix<double, 6, 6>;
using Mat62 = Eigen::Matrix<double, 6, 2>;
using Mat65 = Eigen::Matrix<double, 6, 5>;

enum class Boundary { Pinned, Free };

struct BendingRodSegmentsParameters {
  int links = 14;
  double length = 1.0;
  double total_mass = 1.0;
  double ei = 0.0;
  Boundary boundary = Boundary::Pinned;
  Vec3 gravity = Vec3(0.0, -9.81, 0.0);
  Vec3 pin = Vec3::Zero();

  std::vector<double> masses;
  std::vector<double> inertias;
  std::vector<double> bending_stiffness;
  std::vector<double> bending_damping;

};

struct BendingRodSegmentsConfig : BendingRodSegmentsParameters {
  // Zero uses the complete coupling timestep; positive values select solid
  // substeps (or the initial step when adaptivity is enabled).
  double integration_step = 0.0;
  bool adaptive = false;
  double relative_tolerance = 1e-6;
  double absolute_tolerance = 1e-9;
};

/** A three-dimensional chain of inextensible segments with bending joints. */
class BendingRodSegments : public IBForceCoupled {
public:
  explicit BendingRodSegments(const BendingRodSegmentsConfig &config);

  int links() const;
  int stateSize() const;
  bool isFree() const;

  // Replace the internally owned state and refresh both IB meshes.
  void setState(const Eigen::VectorXd &state);
  const Eigen::VectorXd &state() const { return state_; }
  // Input is force per unit length at the links()+1 endpoints, held constant
  // over dt. Advances the owned state and refreshes GetCurrent()/GetNext().
  void ComputeNextPoints(std::vector<IBMesh::IBVertex> force, real_t dt) override;
  void pack_state(IBModelState &s) const override;
  void unpack_state(const IBModelState &s) override;

  Eigen::VectorXd initialState(const std::vector<Vec3> &tangents = {},
                              const std::vector<Vec3> &angular_velocity = {},
                              Vec3 center_of_mass = Vec3::Zero(),
                              Vec3 center_velocity = Vec3::Zero()) const;
  // Optional endpoint forces are integrated forces, not force densities.
  Eigen::VectorXd rhs(double time, const Eigen::VectorXd &state,
                      const std::vector<Vec3> &nodal_forces = {}) const;
  Eigen::MatrixXd nodes(const Eigen::VectorXd &state) const;
  double energy(const Eigen::VectorXd &state) const;
  void normalizeState(Eigen::VectorXd &state) const;

protected:
  virtual Eigen::VectorXd integrateState(const Eigen::VectorXd &state, double dt,
                                         const std::vector<Vec3> &loads) const;

  struct StateParts {
    Vec3 center = Vec3::Zero();
    Vec3 center_velocity = Vec3::Zero();
    std::vector<Vec3> tangents;
    std::vector<Vec3> omega;
  };

  static void validateNonnegative(double value, const char *name);
  std::vector<Vec3> broadcastLinkValues(const std::vector<Vec3> &supplied,
                                      const Vec3 &fallback,
                                      const char *name) const;
  int tangentOffset() const;
  int velocityOffset() const;
  int omegaOffset() const;
  StateParts unpack(const Eigen::VectorXd &state) const;
  std::vector<Vec3> jointTorques(const std::vector<Vec3> &tangents,
                               const std::vector<Vec3> &omega) const;
  std::vector<Vec3> angularAcceleration(const std::vector<Vec3> &tangents,
                                      const std::vector<Vec3> &omega,
                                      const std::vector<Vec3> &nodal_forces) const;
  void syncMeshes();

  int n_;
  double length_;
  double h_ = 0.0;
  double total_mass_;
  Boundary boundary_;
  Vec3 gravity_;
  Vec3 pin_;
  std::vector<double> masses_;
  std::vector<double> inertias_;
  std::vector<double> bending_stiffness_;
  std::vector<double> bending_damping_;
  Eigen::VectorXd com_weights_;
private:
  Eigen::VectorXd state_;
  double integration_step_;
  bool adaptive_;
  double relative_tolerance_;
  double absolute_tolerance_;
};

// Dormand--Prince fifth-order steps. By default initial_step is fixed, with
// a shorter final step to reach duration. Set adaptive=true to enable the
// embedded fourth-order error estimate and tolerance-based step selection.
Eigen::VectorXd integrateRK45(const BendingRodSegments &beam, Eigen::VectorXd state,
                             double duration, double initial_step,
                             double relative_tolerance = 1e-6,
                             double absolute_tolerance = 1e-9,
                             bool adaptive = false,
                             const std::vector<Vec3> &nodal_forces = {});

} // namespace Models
} // namespace ELFF
