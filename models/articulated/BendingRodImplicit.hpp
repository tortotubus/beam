#pragma once

#include "elff/models/articulated/BendingRod.hpp"
#include <Eigen/Dense>
#include "elff/models/runtime/Model.hpp"

namespace ELFF {
namespace Models {

using Vec3 = Eigen::Vector3d;
using Boundary = BendingRod::Boundary;

/**
 * @brief Geometry, material properties, boundary conditions and integration controls.
 *
 * The rod has \f(n\f) equal-length links and \f(n+1\f) nodes, with
 * \f(h=L/n\f). Empty material arrays select uniform properties; supplied
 * arrays must contain one entry per link or per interior joint, as indicated.
 * All quantities use a consistent system of units.
 */
struct BendingRodImplicitConfig {
  /**
   * @brief Number of rigid links \f(n\geq 1\f); the mesh contains \f(n+1\f) nodes.
   */
  int links = 14;
  /**
   * @brief Total rod length \f(L>0\f).
   */
  real_t length = 1.0;
  /**
   * @brief Positive total mass \f(M\f) used when @ref masses is empty.
   */
  real_t total_mass = 1.0;
  /**
   * @brief Nonnegative bending rigidity \f(EI\f), used for default joint stiffness.
   */
  real_t ei = 0.0;
  /**
   * @brief Fix node zero for Pinned; allow center-of-mass translation for Free.
   */
  Boundary boundary = Boundary::Pinned;
  /**
   * @brief Fixed position of node zero for a pinned rod; unused for a free rod.
   */
  Vec3 pin = Vec3::Zero();

  /**
   * @brief Positive link masses \f(m_i\f); defaults to \f(M/n\f) on each link.
   */
  std::vector<real_t> masses;
  /**
   * @brief Positive transverse inertias \f(I_i\f); defaults to \f(m_i h^2/12\f).
   */
  std::vector<real_t> inertias;
  /**
   * @brief Nonnegative stiffnesses of the \f(n-1\f) joints; defaults to \f(EI/h\f).
   */
  std::vector<real_t> bending_stiffness;
  /**
   * @brief Uniform bending damping rigidity \f(D\geq0\f), in force-length squared-time units.
   * Default joint damping is \f(D/h\f); zero disables physical damping.
   * Explicit @ref bending_damping coefficients override this uniform value.
   */
  real_t damping_rigidity = 0.0;
  /**
   * @brief Nonnegative rotational damping coefficients of the \f(n-1\f) joints; defaults to @ref damping_rigidity divided by link length.
   */
  std::vector<real_t> bending_damping;

  /**
   * @brief High-frequency spectral radius \f(\rho_\infty\in[0,1]\f).
   * Smaller values increase numerical damping of fast linear modes.
   */
  real_t spectral_radius = 0.8;
  /**
   * @brief Positive nonlinear residual tolerance \f(\varepsilon\f).
   * The mass-scaled momentum residual is accepted when
   * \f(\|\mathbf{r}\|_\infty\leq\varepsilon\Delta t
   * (1+\|\mathbf{v}\|_\infty)\f), where \f(\mathbf{v}\f) is the
   * initial velocity in the step's local coordinates.
   */
  real_t nonlinear_tolerance = 1e-9;
  /**
   * @brief Positive maximum number of Newton updates per structural step.
   */
  int max_iterations = 20;
};

/**
 * @brief Implicit dynamics of a spatial, inextensible rod represented by rigid links.
 *
 * Each of the \f(n\f) links has length \f(h=L/n\f), unit tangent
 * \f(\mathbf{t}_i\f), and transverse angular velocity
 * \f(\boldsymbol{\omega}_i\f). Nodal positions satisfy
 * \f[
 *   \mathbf{x}_{i+1}=\mathbf{x}_i+h\mathbf{t}_i,\qquad
 *   \|\mathbf{t}_i\|^2=1,\qquad
 *   \dot{\mathbf{t}}_i=\boldsymbol{\omega}_i\times\mathbf{t}_i,
 *   \qquad \boldsymbol{\omega}_i\cdot\mathbf{t}_i=0.
 * \f]
 * Interior joints resist changes of angle through the bending energy
 * \f[
 *   U_b=\frac{1}{2}\sum_{i=0}^{n-2} k_i\theta_i^2,\qquad
 *   \theta_i=\operatorname{atan2}
 *   \left(\|\mathbf{t}_i\times\mathbf{t}_{i+1}\|,
 *   \mathbf{t}_i\cdot\mathbf{t}_{i+1}\right).
 * \f]
 * The model includes transverse rotational inertia and optional joint damping,
 * but no axial extension or torsional degree of freedom. A pinned boundary
 * fixes node zero while leaving the first link free to rotate. A free boundary
 * adds center-of-mass position and velocity to the state.
 *
 * Time integration uses a generalized-alpha step in local two-dimensional
 * coordinates on the unit sphere for each tangent. With spectral radius
 * \f(\rho_\infty\f), the implementation uses
 * \f[
 *   \alpha_m=\frac{2-\rho_\infty}{1+\rho_\infty},\qquad
 *   \alpha_f=\frac{1}{1+\rho_\infty},\qquad
 *   \gamma=\frac{1}{2}+\alpha_m-\alpha_f,\qquad
 *   \beta=\frac{1}{4}(1+\alpha_m-\alpha_f)^2.
 * \f]
 * Here \f(\alpha_m\f) and \f(\alpha_f\f) weight the new acceleration
 * and the new configuration/velocity, respectively. External loads remain
 * fixed during each step. Numerical acceleration history is reused only when
 * the timestep matches the previous step.
 *
 * @note Public solve() loads are force densities per unit length. Protected
 * advance() and rhs() take integrated nodal forces.
 * @note Invalid inputs and nonlinear failure use ELFF_ABORT and follow the
 * configured ELFF error policy. Rejected steps do not commit a new state or
 * acceleration history.
 */
class BendingRodImplicit : public BendingRod {
public:
  /**
   * @brief Construct a straight rod at rest and initialize the material operators.
   * @param config Geometry, material, boundary and generalized-alpha controls.
   * The initial tangents point along the positive x-axis. Node zero is at the
   * configured pin for a pinned rod; the center of mass is at the origin for a
   * free rod. Use apply_initial_condition() to prescribe another initial mesh.
   */
  explicit BendingRodImplicit(const BendingRodImplicitConfig &config);
  using BendingRod::solve;
  /**
   * @brief Advance one structural step with prescribed nodal force densities.
   * @param dt Finite, nonnegative timestep \f(\Delta t\f); zero leaves the state unchanged.
   * @param load One finite force density \f(\mathbf{q}_j\f) per mesh node.
   * The integrated force is \f(\mathbf{F}_j=w_j\mathbf{q}_j\f), with
   * \f(w_0=w_n=h/2\f) and \f(w_j=h\f) at interior nodes. A successful
   * step updates the common mesh positions and velocities.
   */
  void solve(real_t dt, const std::vector<Vector>& load) override;
  /**
   * @brief Initialize the structural state from nodal positions and velocities.
   * @param mesh Mesh with the same node count and length as the solver.
   * The mesh must have finite data, link lengths \f(h\f), and velocities
   * satisfying \f(\mathbf{t}_i\cdot(\mathbf{v}_{i+1}-\mathbf{v}_i)=0\f).
   * For a pinned rod, node zero must remain at the pin with zero velocity.
   * Initialization resets the acceleration history and previous timestep.
   */
  void apply_initial_condition(const BendingRodMesh& mesh) override;
  /**
   * @brief Serialize structural parameters, state and generalized-alpha history.
   * @param[out] s Checkpoint container whose previous contents are replaced.
   * The structural payload excludes coupling forces; an IB adapter may append them.
   */
  virtual void pack_state(ModelState &s) const;
  /**
   * @brief Restore a compatible checkpoint and refresh the common mesh.
   * @param s Structural checkpoint produced by pack_state().
   * Model parameters and payload dimensions must match this solver. Restoration
   * preserves the saved acceleration history and timestep rather than resetting them.
   */
  virtual void unpack_state(const ModelState &s);
  /**
   * @brief Return the current kinetic plus bending energy, excluding external-load potential.
   * @return Mechanical energy; excludes external-load work and dissipated energy.
   */
  real_t energy() const { return energy(state_); }

protected:
  /**
   * @brief Reconstruct the common mesh positions and velocities from @ref state_.
   */
  void update_mesh();
  /**
   * @brief Compute one generalized-alpha step using a Newton solve with line search.
   * @param state Structural state at the beginning of the step.
   * @param dt Finite, strictly positive timestep.
   * @param loads Integrated nodal forces, held fixed throughout the step; empty means zero.
   * @return Updated structural state, without assigning @ref state_ or updating the mesh.
   * On success this method updates the mutable acceleration history and previous
   * timestep. On validation or convergence failure, that history is unchanged.
   * Logs one convergence summary through ELFF_LOG: accepted Newton updates,
   * residual evaluations, rejected line-search trials, and initial/final residual
   * infinity norms against the timestep-scaled stopping tolerance. Zero updates
   * means the initial guess already satisfies that tolerance.
   */
  Eigen::VectorXd integrateState(const Eigen::VectorXd &state, real_t dt,
                                const std::vector<Vec3> &loads) const;
  /**
   * @brief Validate integrated nodal loads, advance the state and refresh the mesh.
   * @param dt Finite, nonnegative timestep; zero validates inputs without advancing.
   * @param nodal_forces Integrated forces on all \f(n+1\f) nodes; empty means zero.
   * The new state is committed after integrateState() succeeds.
   */
  void advance(real_t dt, const std::vector<Vec3>& nodal_forces = {});
  /**
   * @brief Return the stored structural state without modifying its constraints or history.
   */
  const Eigen::VectorXd &state() const { return state_; }
  /**
   * @brief Validate and normalize a replacement state, then refresh the mesh.
   * @param state State vector in the layout documented by @ref state_.
   * Clears numerical acceleration history. Derived adapters may also synchronize
   * their coupling meshes and reset coupling forces.
   */
  virtual void setState(const Eigen::VectorXd &state);
  /**
   * @brief Return \f(6n\f) for a pinned rod or \f(6n+6\f) for a free rod.
   */
  int stateSize() const;
  /**
   * @brief Return whether center-of-mass translation is an independent degree of freedom.
   */
  bool isFree() const;
  /**
   * @brief Reconstruct nodal velocities from a structural state.
   * @param state Structural state in the layout documented by @ref state_.
   * @return An \f((n+1)\times3\f) matrix, one Cartesian velocity per row.
   * Adjacent velocities satisfy \f(\mathbf{v}_{i+1}-\mathbf{v}_i
   * =h\boldsymbol{\omega}_i\times\mathbf{t}_i\f). The translational
   * shift enforces the fixed pin or the prescribed center-of-mass velocity.
   */
  Eigen::MatrixXd velocities(const Eigen::VectorXd &state) const;
  /**
   * @brief Construct a normalized structural state from link data.
   * @param tangents Empty for the positive x-axis, one vector to broadcast, or \f(n\f) vectors.
   * @param angular_velocity Empty for rest, one vector to broadcast, or \f(n\f) vectors.
   * @param center_of_mass Center-of-mass position for a free rod; ignored when pinned.
   * @param center_velocity Center-of-mass velocity for a free rod; ignored when pinned.
   * @return State with unit tangents and angular velocities projected perpendicular to them.
   * This method constructs a state without assigning it or resetting solver history.
   */
  Eigen::VectorXd initialState(const std::vector<Vec3> &tangents = {},
                              const std::vector<Vec3> &angular_velocity = {},
                              Vec3 center_of_mass = Vec3::Zero(),
                              Vec3 center_velocity = Vec3::Zero()) const;
  /**
   * @brief Evaluate the continuous structural equations of motion.
   * @param time Time argument, unused because the supplied loads are instantaneous.
   * @param state Structural state at which to evaluate the derivative.
   * @param nodal_forces Integrated nodal forces; empty means zero external load.
   * @return State derivative including applied nodal forces, bending and damping. For a free rod,
   * \f(\ddot{\mathbf{c}}=M^{-1}\sum_j\mathbf{F}_j\f).
   */
  Eigen::VectorXd rhs(real_t time, const Eigen::VectorXd &state,
                      const std::vector<Vec3> &nodal_forces = {}) const;
  /**
   * @brief Reconstruct node positions from the link tangents and boundary condition.
   * @param state Structural state in the layout documented by @ref state_.
   * @return An \f((n+1)\times3\f) matrix, one Cartesian position per row.
   * The first node is the pin or is inferred from the free center of mass.
   */
  Eigen::MatrixXd nodes(const Eigen::VectorXd &state) const;
  /**
   * @brief Evaluate mechanical energy for a supplied structural state.
   * @param state Structural state at which to evaluate energy.
   * @return Kinetic plus bending energy \f(T+U_b\f).
   * External-load potential energy is not included.
   */
  real_t energy(const Eigen::VectorXd &state) const;
  /**
   * @brief Enforce unit tangents and remove axial spin from a state in place.
   * @param[in,out] state Vector with the expected state size and nonzero finite tangents.
   */
  void normalizeState(Eigen::VectorXd &state) const;

  /**
   * @brief Unpacked kinematic fields with normalized tangents and transverse angular velocities.
   */
  struct StateParts {
    /**
     * @brief Center-of-mass position \f(\mathbf{c}\f) for a free rod.
     */
    Vec3 center = Vec3::Zero();
    /**
     * @brief Center-of-mass velocity \f(\dot{\mathbf{c}}\f) for a free rod.
     */
    Vec3 center_velocity = Vec3::Zero();
    /**
     * @brief Unit link tangents \f(\mathbf{t}_i\f).
     */
    std::vector<Vec3> tangents;
    /**
     * @brief Transverse link angular velocities \f(\boldsymbol{\omega}_i\f).
     */
    std::vector<Vec3> omega;
  };

  /**
   * @brief Require a finite, nonnegative value, using the supplied name in diagnostics.
   */
  static void validateNonnegative(real_t value, const char *name);
  /**
   * @brief Expand an empty array to the fallback, broadcast one value, or accept exactly \f(n\f) values.
   */
  std::vector<Vec3> broadcastLinkValues(const std::vector<Vec3> &supplied,
                                      const Vec3 &fallback,
                                      const char *name) const;
  /**
   * @brief Return the first tangent component index in the structural state.
   */
  int tangentOffset() const;
  /**
   * @brief Return the first velocity component index, after all tangent components.
   */
  int velocityOffset() const;
  /**
   * @brief Return the first angular velocity component index, after any center velocity.
   */
  int omegaOffset() const;
  /**
   * @brief Validate a state and extract normalized kinematics without modifying the input.
   */
  StateParts unpack(const Eigen::VectorXd &state) const;
  /**
   * @brief Assemble equal-and-opposite elastic and damping torques at interior joints.
   * @param tangents Unit link tangents.
   * @param omega Transverse link angular velocities.
   * @return Net joint torque on each link; damping acts on relative angular
   * velocity projected perpendicular to the joint mean tangent.
   */
  std::vector<Vec3> jointTorques(const std::vector<Vec3> &tangents,
                               const std::vector<Vec3> &omega) const;
  /**
   * @brief Compute link angular accelerations by articulated-body elimination.
   * @param tangents Unit link tangents.
   * @param omega Transverse link angular velocities.
   * @param nodal_forces Integrated nodal forces; empty means zero external load.
   * @return Angular acceleration of each link, accounting for nodal forces, joint
   * forces and centripetal terms. The pinned base has two rotational freedoms;
   * the free base additionally has three translational freedoms.
   */
  std::vector<Vec3> angularAcceleration(const std::vector<Vec3> &tangents,
                                      const std::vector<Vec3> &omega,
                                      const std::vector<Vec3> &nodal_forces) const;

  /**
   * @brief Number of links \f(n\f).
   */
  int n_;
  /**
   * @brief Total length \f(L\f).
   */
  real_t length_;
  /**
   * @brief Uniform link length \f(h=L/n\f).
   */
  real_t h_ = 0.0;
  /**
   * @brief Sum of the resolved link masses \f(M=\sum_i m_i\f).
   */
  real_t total_mass_;
  /**
   * @brief Pinned or free translational boundary condition.
   */
  Boundary boundary_;
  /**
   * @brief Fixed first-node position for a pinned rod.
   */
  Vec3 pin_;
  /**
   * @brief Resolved positive mass \f(m_i\f) of each link.
   */
  std::vector<real_t> masses_;
  /**
   * @brief Resolved positive transverse inertia \f(I_i\f) of each link.
   */
  std::vector<real_t> inertias_;
  /**
   * @brief Resolved nonnegative stiffness \f(k_i\f) of each interior joint.
   */
  std::vector<real_t> bending_stiffness_;
  /**
   * @brief Resolved nonnegative damping coefficient of each interior joint.
   */
  std::vector<real_t> bending_damping_;
  /**
   * @brief Tangent weights locating the center of mass relative to node zero.
   * \f[ c_i=\frac{h}{M}\left(\frac{m_i}{2}+\sum_{j=i+1}^{n-1}m_j\right),
   * \qquad \mathbf{c}=\mathbf{x}_0+\sum_i c_i\mathbf{t}_i. \f]
   */
  Eigen::VectorXd com_weights_;

  /**
   * @brief Packed structural state in Cartesian components.
   * Pinned layout: \f([\mathbf{t}_0,\ldots,\mathbf{t}_{n-1},
   * \boldsymbol{\omega}_0,\ldots,\boldsymbol{\omega}_{n-1}]\f).
   * Free layout: \f([\mathbf{c},\mathbf{t}_0,\ldots,\mathbf{t}_{n-1},
   * \dot{\mathbf{c}},\boldsymbol{\omega}_0,\ldots,\boldsymbol{\omega}_{n-1}]\f).
   */
  Eigen::VectorXd state_;
  /**
   * @brief Constant \f(n\times n\f) mass matrix for ambient tangent coordinates.
   * Includes link translational mass and transverse inertia. For a free rod,
   * the translational contribution is expressed relative to the center of mass.
   */
  Eigen::MatrixXd tangent_mass_;
  /**
   * @brief Generalized-alpha spectral radius \f(\rho_\infty\f) and nonlinear
   * residual tolerance \f(\varepsilon\f), respectively.
   */
  real_t spectral_radius_, tolerance_;
  /**
   * @brief Maximum number of Newton updates per structural step.
   */
  int max_iterations_;
  /**
   * @brief Numerical acceleration history in ambient coordinates.
   * Stores \f(\ddot{\mathbf{t}}_i\f) for each link, preceded by
   * \f(\ddot{\mathbf{c}}\f) for a free rod. Empty until a step succeeds
   * or after explicit reinitialization; updated by integrateState().
   */
  mutable Eigen::VectorXd acceleration_;
  /**
   * @brief Timestep associated with @ref acceleration_; zero when history is absent.
   * A changed timestep causes the next step to initialize acceleration from rhs().
   */
  mutable real_t previous_dt_ = 0.0;
};

} // namespace Models
} // namespace ELFF
