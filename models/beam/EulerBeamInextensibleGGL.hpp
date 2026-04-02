// ============================================================
// EulerBeamInextensibleGGL.hpp
// ============================================================
#pragma once

#include "elff/models/beam/EulerBeam.hpp"
#include "elff/fem/Shapes.hpp"

#include <array>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <unsupported/Eigen/AutoDiff>

using namespace Eigen;

namespace ELFF {
namespace Models {

/**
 * @brief Dynamic inextensible Euler beam solved with a full
 * Gear-Gupta-Leimkuhler system in the independent unknowns
 * (u, v, lambda, mu).
 *
 * The nonlinear unknown at each Newton step is stacked as
 *
 *   [ u (ndof) | v (ndof) | lambda (ndof_l) | mu (ndof_l) ].
 *
 * The four residual blocks are:
 *
 *   R_u      = M a(u) + K(u, lambda) - f
 *   R_v      = M (v - v_newmark(u)) - C(u)^T mu
 *   R_lambda = g(u)
 *   R_mu     = C(u) v
 *
 * where:
 *   a(u)          = 1/(beta dt^2) (u - u_tilde)
 *   v_newmark(u)  = v_tilde + gamma/(beta dt) (u - u_tilde)
 *   g(u)          = ||r'||^2 - 1
 *   C(u)          = d/ dv [ G(u) v ]
 *
 * This keeps the velocity constraint in its own variable space,
 * so no Schur-complement mass inversion is required.
 */
class EulerBeamInextensibleGGL
  : public EulerBeam
{
public:
  /**
   * @brief Construct a dynamic inextensible Euler beam with a full GGL solve.
   *
   * @param length Total beam length.
   * @param EI Bending stiffness.
   * @param mu_density Mass per unit length.
   * @param n_nodes Number of mesh nodes.
   * @param bcs Boundary conditions inherited from the base beam model.
   * @param r_penalty Optional static regularization penalty passed to the base
   *   class.
   */
  EulerBeamInextensibleGGL(real_t                  length,
                                  real_t                  EI,
                                  real_t                  mu_density,
                                  size_t                  n_nodes,
                                  EulerBeam::EulerBeamBCs bcs,
                                  real_t                  r_penalty = 0.0);

  ~EulerBeamInextensibleGGL() = default;

  /**
   * @brief Advance one time step under a spatially uniform load.
   *
   * @param dt Time-step size.
   * @param load Uniform body force per unit length.
   */
  void solve(real_t dt, std::array<real_t, 3> load) override;
  /**
   * @brief Advance one time step under nodal load data.
   *
   * The nodal values are linearly interpolated elementwise during assembly.
   *
   * @param dt Time-step size.
   * @param load Nodal load values with one vector per mesh node.
   */
  void solve(real_t dt, std::vector<std::array<real_t, 3>> load) override;

  /**
   * @brief Reset the dynamic state from the current beam configuration.
   */
  void apply_initial_condition();
  /**
   * @brief Reset the dynamic state from an externally supplied beam mesh.
   *
   * The mesh contributes both the centerline geometry and the initial
   * centerline velocity.
   *
   * @param bmesh Mesh used to seed the dynamic state.
   */
  void apply_initial_condition(EulerBeamMesh& bmesh) override;

  /**
   * @brief Return the converged physical beam tension field.
   *
   * In the present formulation the physical tension is represented by the
   * converged position multiplier.
   */
  const VectorXd& get_tension() const { return tension; }

  /**
   * @brief Return the current physical time.
   */
  real_t get_time() const { return t; }
  /**
   * @brief Return the completed time-step count.
   */
  size_t get_time_iter() const { return time_iter; }

protected:
  size_t elements; ///< Number of two-node Hermite beam elements.
  size_t nodes;    ///< Number of mesh nodes.
  real_t ds;       ///< Uniform nodal spacing along the centerline.
  size_t ndof_l;   ///< Number of multiplier degrees of freedom.
  size_t offset_x; ///< Global offset of the x block in displacement-like vectors.
  size_t offset_y; ///< Global offset of the y block in displacement-like vectors.
  size_t offset_z; ///< Global offset of the z block in displacement-like vectors.
  size_t ndof;     ///< Total number of displacement-like degrees of freedom.
  real_t r_penalty; ///< Retained penalty parameter for interface compatibility.

  VectorXd u;      ///< Current displacement/slope state packed by component.
  VectorXd u_prev;      ///< Displacement state at the previous time step.
  VectorXd v_prev;      ///< Velocity state at the previous time step.
  VectorXd a_prev;      ///< Acceleration state at the previous time step.
  VectorXd lambda;      ///< Current position multiplier field.
  VectorXd lambda_prev; ///< Previous-step position multiplier.
  VectorXd mu_prev;     ///< Previous-step velocity multiplier.
  VectorXd tension;     ///< Converged physical tension field for output.

  SparseMatrix<real_t> mass; ///< Consistent Hermite mass matrix.

  size_t offset_v;      ///< Global offset of the velocity block.
  size_t offset_lambda; ///< Global offset of the position multiplier block.
  size_t offset_mu;     ///< Global offset of the velocity multiplier block.
  size_t ggl_dof;       ///< Total nonlinear system size.

  VectorXd             ggl_residual; ///< Full stacked GGL residual.
  SparseMatrix<real_t> ggl_jacobian; ///< Full stacked GGL Jacobian.

  real_t newmark_coeff_a; ///< Newmark coefficient mapping displacement to acceleration.
  real_t newmark_coeff_v; ///< Newmark coefficient mapping displacement to velocity.
  real_t newmark_gamma;   ///< Newmark gamma parameter.
  real_t newmark_beta;    ///< Newmark beta parameter.
  VectorXd u_tilde;       ///< Newmark displacement predictor.
  VectorXd v_tilde;       ///< Newmark velocity predictor.

  std::array<std::array<real_t, 4>, 3> quad_H;   ///< Cached Hermite values at quadrature points.
  std::array<std::array<real_t, 4>, 3> quad_dH;  ///< Cached Hermite first derivatives.
  std::array<std::array<real_t, 4>, 3> quad_ddH; ///< Cached Hermite second derivatives.
  std::array<std::array<real_t, 2>, 3> quad_M;   ///< Cached linear multiplier shape values.
  std::vector<std::array<size_t, 12>>
    element_disp_dof_indices_cache; ///< Cached displacement dof maps per element.
  Matrix<real_t, 12, 12>
    local_bending_jacobian; ///< Constant local bending tangent reused during assembly.

  size_t max_newton; ///< Maximum Newton iterations per time step.
  real_t tol_newton; ///< Newton stopping tolerance on the full residual norm.

  /** @brief Precompute all shape-function data used at the fixed quadrature rule. */
  void initialize_quadrature_cache();
  /** @brief Cache the global displacement dof indices for each element. */
  void initialize_element_dof_cache();
  /** @brief Build constant element matrices used by every Newton iteration. */
  void initialize_constant_element_matrices();
  /** @brief Assemble the consistent translational-and-slope mass matrix. */
  void assemble_mass_matrix();

  /**
   * @brief Assemble the full GGL residual and Jacobian for a uniform load.
   *
   * @param u_cur Current displacement iterate.
   * @param v_cur Current velocity iterate.
   * @param lambda_cur Current position multiplier iterate.
   * @param mu_cur Current velocity multiplier iterate.
   * @param load Uniform body force per unit length.
   */
  void assemble_ggl_system(const VectorXd& u_cur,
                           const VectorXd& v_cur,
                           const VectorXd& lambda_cur,
                           const VectorXd& mu_cur,
                           std::array<real_t, 3> load);

  /**
   * @brief Assemble the full GGL residual and Jacobian for nodal load data.
   *
   * @param u_cur Current displacement iterate.
   * @param v_cur Current velocity iterate.
   * @param lambda_cur Current position multiplier iterate.
   * @param mu_cur Current velocity multiplier iterate.
   * @param load Nodal load data interpolated elementwise during assembly.
   */
  void assemble_ggl_system(const VectorXd&                            u_cur,
                           const VectorXd&                            v_cur,
                           const VectorXd&                            lambda_cur,
                           const VectorXd&                            mu_cur,
                           const std::vector<std::array<real_t, 3>>& load);

  /**
   * @brief Assemble one element residual for a uniform distributed load.
   *
   * The returned 16-vector is ordered as [R_u(12), R_lambda(2), R_mu(2)].
   */
  template<typename Scalar>
  Matrix<Scalar, 16, 1> assemble_ggl_element(
    const Matrix<Scalar, 12, 1>& u_elem,
    const std::array<real_t, 2>& lambda_elem,
    const std::array<real_t, 4>& vx_elem,
    const std::array<real_t, 4>& vy_elem,
    const std::array<real_t, 4>& vz_elem,
    std::array<real_t, 3> load) const;

  /**
   * @brief Assemble one element residual for linearly interpolated nodal load
   * data.
   *
   * The returned 16-vector is ordered as [R_u(12), R_lambda(2), R_mu(2)].
   */
  template<typename Scalar>
  Matrix<Scalar, 16, 1> assemble_ggl_element(
    const Matrix<Scalar, 12, 1>&                u_elem,
    const std::array<real_t, 2>&                lambda_elem,
    const std::array<real_t, 4>&                vx_elem,
    const std::array<real_t, 4>&                vy_elem,
    const std::array<real_t, 4>&                vz_elem,
    const std::array<std::array<real_t, 3>, 2>& load_elem) const;

  /**
   * @brief Assemble the velocity-row contribution induced by the multiplier mu.
   *
   * This is the local term -C(u)^T mu appearing in the discrete velocity
   * equation.
   */
  template<typename Scalar>
  Matrix<Scalar, 12, 1> assemble_velocity_multiplier_residual(
    const Matrix<Scalar, 12, 1>& u_elem,
    const std::array<real_t, 2>& mu_elem) const;

  /** @brief Collect constrained dofs and their prescribed values. */
  void collect_ggl_constraints(std::vector<char>& constrained,
                               VectorXd&          constrained_value) const;

  /** @brief Update Newmark predictor states and coefficients for the next solve. */
  void update_newmark_predictors(real_t dt, real_t beta, real_t gamma);

  /** @brief Recover the Newmark-consistent velocity for a displacement iterate. */
  VectorXd compute_velocity(const VectorXd& u_cur) const;
  /** @brief Recover the Newmark-consistent acceleration for a displacement iterate. */
  VectorXd compute_acceleration(const VectorXd& u_cur) const;

  /** @brief Compute quadrature-based inextensibility errors for a state vector. */
  std::pair<real_t, real_t> compute_inextensibility_error(
    const VectorXd& u_cur) const;

  /** @brief Commit the converged displacement and velocity state. */
  void update_newmark_state(const VectorXd& u_new, const VectorXd& v_new);

  /** @brief Reapply essential boundary conditions to the stored dynamic state. */
  void apply_dynamic_state_boundary_conditions();

  /** @brief Push the converged solution back into the beam mesh representation. */
  void update_mesh();

  /** @brief Return the 12 local displacement DOF indices for one element. */
  const std::array<size_t, 12>& get_element_disp_dof_indices(size_t e) const;

  /** @brief Extract the element x-velocity DOFs from a global vector. */
  std::array<real_t, 4> get_element_velocity_x(size_t e,
                                               const VectorXd& v) const;
  /** @brief Extract the element y-velocity DOFs from a global vector. */
  std::array<real_t, 4> get_element_velocity_y(size_t e,
                                               const VectorXd& v) const;
  /** @brief Extract the element z-velocity DOFs from a global vector. */
  std::array<real_t, 4> get_element_velocity_z(size_t e,
                                               const VectorXd& v) const;
};

} // namespace Models
} // namespace ELFF
