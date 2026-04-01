#pragma once

#include "elff/models/beam/EulerBeamInextensibleMoM.hpp"

using namespace Eigen;

namespace ELFF {
namespace Models {

/**
 * @brief Dynamic inextensible Euler beam solved as an index-1 DAE with a
 * sparse saddle-point system and Newmark time integration.
 *
 * The inextensibility constraint is enforced at the velocity level, yielding
 * an index-1 system.  Drift in the position and velocity constraints is
 * corrected via projection after each time step.
 */
class EulerBeamDynamicInextensibleIndex1 : public EulerBeamInextensibleMoM
{
public:
  /**
   * @brief Constructs the index-1 dynamic inextensible sparse beam model.
   *
   * @param length   Beam length
   * @param EI       Flexural rigidity
   * @param mu       Mass per unit length
   * @param n_nodes  Number of discretization nodes
   * @param bcs      Boundary conditions at the beam ends
   */
  EulerBeamDynamicInextensibleIndex1(real_t                  length,
                                     real_t                  EI,
                                     real_t                  mu,
                                     size_t                  n_nodes,
                                     EulerBeam::EulerBeamBCs bcs);

  /**
   * @brief Applies the model's default initial condition.
   */
  void apply_initial_condition() override;

  /**
   * @brief Applies initial conditions from a supplied beam mesh.
   *
   * @param bmesh Beam mesh containing the initial geometry and velocities
   */
  void apply_initial_condition(EulerBeamMesh& bmesh) override;

  /**
   * @brief Advances the beam one time step under a uniform load.
   *
   * @param dt   Time-step size
   * @param load Uniform distributed load vector
   */
  void solve(real_t dt, std::array<real_t, 3> load) override;

  /**
   * @brief Advances the beam one time step under nodal loads.
   *
   * @param dt   Time-step size
   * @param load Load vector specified at the mesh nodes
   */
  void solve(real_t dt, std::vector<std::array<real_t, 3>> load) override;

protected:
  /** @brief Displacement state from the previous time step. */
  VectorXd u_prev;
  /** @brief Velocity degrees of freedom from the previous time step. */
  VectorXd v_prev;
  /** @brief Acceleration degrees of freedom from the previous time step. */
  VectorXd a_prev;

  /** @brief Constant elastic stiffness matrix (assembled once). */
  SparseMatrix<real_t> K_elastic;
  /** @brief Saddle-point system matrix assembled each time step. */
  SparseMatrix<real_t> saddle_mat;
  /** @brief Right-hand side for the saddle-point system. */
  VectorXd saddle_rhs;

  /** @brief Tolerance for position constraint drift correction. */
  real_t tol_position_drift;
  /** @brief Tolerance for velocity constraint drift correction. */
  real_t tol_velocity_drift;
  /** @brief Maximum number of position projection iterations. */
  size_t max_projection_iter;

  /** @brief Current time-step index. */
  size_t time_iter;
  /** @brief Current simulation time. */
  real_t t;

  /**
   * @brief Assembles the constant elastic stiffness matrix K_elastic.
   */
  void assemble_elastic_stiffness();

  /**
   * @brief Assembles the constraint Jacobian B and the velocity-level RHS g.
   *
   * B has size ndof_l x ndof; g has size ndof_l.
   *
   * @param B Output constraint Jacobian
   * @param g Output velocity-level right-hand side
   */
  void assemble_B_and_g(SparseMatrix<real_t>& B, VectorXd& g) const;

  /**
   * @brief Assembles the position constraint residual vector.
   *
   * @return Vector c where c_a = int eta_a * (||r'||^2 - 1) ds
   */
  VectorXd assemble_constraint_residual_vector() const;

  /**
   * @brief Assembles the external force vector for a uniform load.
   *
   * @param load Uniform distributed load vector
   * @return Global external force vector
   */
  VectorXd assemble_f_ext(std::array<real_t, 3> load) const;

  /**
   * @brief Assembles the external force vector for a nodal load.
   *
   * @param load Load vector specified at the mesh nodes
   * @return Global external force vector
   */
  VectorXd assemble_f_ext(const std::vector<std::array<real_t, 3>>& load) const;

  /**
   * @brief Assembles the saddle-point system matrix and right-hand side.
   *
   * @param dt    Time-step size
   * @param beta  Newmark beta parameter
   * @param f_ext External force vector
   * @param B     Constraint Jacobian
   * @param g     Velocity-level right-hand side
   */
  void assemble_saddle_system(real_t                      dt,
                               real_t                      beta,
                               const VectorXd&             f_ext,
                               const SparseMatrix<real_t>& B,
                               const VectorXd&             g);

  /**
   * @brief Enforces displacement and tension boundary conditions on the
   * assembled saddle-point system.
   */
  void apply_saddle_boundary_conditions();

  /**
   * @brief Updates the Newmark velocity and acceleration from the converged
   * displacement.
   *
   * @param dt    Time-step size
   * @param beta  Newmark beta parameter
   * @param gamma Newmark gamma parameter
   */
  void update_newmark_state(real_t dt, real_t beta, real_t gamma);

  /**
   * @brief Zeros velocity and acceleration at kinematically constrained DOFs.
   */
  void apply_dynamic_state_boundary_conditions();

  /**
   * @brief Projects the displacement onto the position constraint manifold
   * via a Newton iteration.
   */
  void project_position_onto_constraint();

  /**
   * @brief Projects the velocity onto the velocity constraint manifold.
   *
   * @param B Constraint Jacobian evaluated at the current configuration
   */
  void project_velocity_onto_constraint(const SparseMatrix<real_t>& B);

  /**
   * @brief Updates the beam mesh positions, slopes, and velocities from the
   * current state vectors.
   */
  void update_mesh();
};

} // namespace Models
} // namespace ELFF
