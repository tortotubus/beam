#pragma once

#include "elff/fem/Shapes.hpp"
#include "elff/models/beam/EulerBeam.hpp"

#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <unsupported/Eigen/AutoDiff>

using namespace Eigen;

namespace ELFF {
namespace Models {

/**
 * @brief A class to solve the inextensible Euler–Bernoulli beam equation in
 * static and time-dependent settings; this model is valid for large
 * deflections. The strong form of the system that we seek to solve is
 * \f[
 *    \frac{\partial^2}{\partial s^2} \left(EI \frac{\partial^2
 * \mathbf{r}}{\partial s^2}\right) = \mathbf{q}(s)
 * \f]
 * where \f(EI\f) is a bending stiffness, \f(\mathbf{r}(s) = (x(s),y(s))\f) is
 * the deflection of our beam. Unlike the classic Euler-Bernouli beam equation,
 * where inextensibility is implicitly enforced, in the \f(n > 1\f) dimensional
 * version, we consider the deflection of the beam in each dimension along a
 * curvilinear coordinate system, and our inextensbility must be enforced
 * explicitly. Then, in addition, we enforce on the solution a pointwise
 * constraint
 * \f[
 *    ||\mathbf{r}'(s)||^2 = 1.
 * \f]
 * To derive a weak form, we introduce a smooth test function
 * ...
 * \f[
 *    \mathcal{L}_R(x,y,p,q,\lambda_x,\mu) = J(x,y) + \int_0^L
 * \left[\lambda_x(p-x')
 * + \mu(q-y')\right] ds + \frac{r}{2} \int_0^L \left[(p-x')^2 +
 * (q-y')^2\right] ds.
 * \f]
 */
class EulerBeamInextensibleADDM : public EulerBeam {
public:
  /**
   * @brief
   */
  EulerBeamInextensibleADDM(real_t length, real_t EI, size_t nodes,
                            EulerBeamBCs bcs, real_t r_penalty);

  /**
   * @brief
   */
  EulerBeamInextensibleADDM(real_t length, real_t EI, real_t mu, size_t nodes,
                            EulerBeamBCs bcs, real_t r_penalty);
  /**
   * @brief
   */
  EulerBeamInextensibleADDM(real_t length, real_t EI, real_t mu, size_t nodes,
                            EulerBeamTimeDependentBCs bcs, real_t r_penalty);

  /**
   * @brief
   */
  ~EulerBeamInextensibleADDM();

  /**
   * @brief
   */
  virtual void solve() override;

  /**
   * @brief
   */
  virtual void solve(std::array<real_t, 3> load) override;

  /**
   * @brief
   */
  virtual void solve(real_t dt) override;

  /**
   * @brief
   */
  virtual void solve(real_t dt, std::array<real_t, 3> load) override;

  /**
   * @brief
   */
  virtual void solve(real_t dt,
                     std::vector<std::array<real_t, 3>> load) override;

  /**
   * @brief Advance one step with a nodal load that is already averaged over
   * the Newmark time interval.
   */
  void solve_averaged_load(
      real_t dt, const std::vector<std::array<real_t, 3>> &averaged_load);

  /**
   * @brief
   */
  void solve_newmark(real_t dt, std::array<real_t, 3> load, real_t beta,
                     real_t gamma);

  /**
   * @brief
   */
  void solve_newmark(real_t dt, std::vector<std::array<real_t, 3>> load,
                     real_t beta, real_t gamma);

  /**
   * @brief
   */
  void solve_newmark_averaged_load(
      real_t dt, const std::vector<std::array<real_t, 3>> &averaged_load,
      real_t beta, real_t gamma);

  /**
   * @brief
   */
  virtual void apply_initial_condition();

  /**
   * @brief
   */
  virtual void apply_initial_condition(EulerBeamMesh &bmesh) override;

  /**
   * @brief
   */
  void set_time_dependent_boundary_conditions(
      EulerBeamTimeDependentBCs bcs) override;

  /**
   * @brief
   */
  void set_outer_tolerance(real_t tol) { tol_outer = tol; }

  /**
   * @brief
   */
  void set_outer_iter_max(int max_iter) { max_outer = max_iter; }

  /**
   * @brief
   */
  void set_outer_relaxation(real_t omega) { omega_outer = omega; }

  /**
   * @brief
   */
  const VectorXd &get_lambda_x() const;

  /**
   * @brief
   */
  const VectorXd &get_lambda_y() const;

  /**
   * @brief
   */
  const VectorXd &get_lambda_z() const;

  /**
   * @brief
   */
  const MatrixXd &get_A() const;

protected:
  size_t dimension;
  size_t elements;
  size_t dof;
  real_t r_penalty, alpha;

  SparseMatrix<real_t> A;
  SparseMatrix<real_t> A_unconstrained;
  SparseMatrix<real_t> K_bending;
  SparseMatrix<real_t> K_constraint;
  mutable MatrixXd A_dense_cache;
  Matrix<real_t, 4, 4> bending_element_matrix;
  Matrix<real_t, 4, 3> constraint_rhs_element_matrix;
  std::array<real_t, 4> midpoint_dH;
  VectorXd x, y, z;
  VectorXd f_x, f_y, f_z;
  SimplicialLLT<SparseMatrix<real_t>> llt;
  // SimplicialLDLT<SparseMatrix<real_t>> llt;

  VectorXd lambda_x, lambda_y, lambda_z;
  VectorXd p, q, r;
  VectorXd xp, yp, zp;

  size_t max_outer, min_outer;
  real_t omega_outer;
  real_t tol_outer;
  VectorXd x_prev;
  VectorXd y_prev;
  VectorXd z_prev;
  VectorXd vx_prev;
  VectorXd vy_prev;
  VectorXd vz_prev;
  VectorXd ax_prev;
  VectorXd ay_prev;
  VectorXd az_prev;
  SparseMatrix<real_t> mass_matrix;
  std::array<real_t, 3> load_prev;
  std::vector<std::array<real_t, 3>> nodal_load_prev;
  bool have_prev_uniform_load;
  bool have_prev_nodal_load;

  /**
   * @brief
   */
  void update_mesh();

  /**
   * @brief
   */
  void apply_time_dependent_boundary_values(size_t step_idx);

  /**
   * @brief
   */
  void apply_time_dependent_boundary_kinematics(size_t step_idx);

  /**
   * @brief
   */
  void apply_prescribed_boundary_values(VectorXd &x_state, VectorXd &y_state,
                                        VectorXd &z_state) const;

  /**
   * @brief
   */
  void clear_constrained_dynamic_history();

  /**
   * @brief
   */
  void collect_boundary_dofs(std::vector<size_t> &idx,
                             std::vector<real_t> &xvals,
                             std::vector<real_t> &yvals,
                             std::vector<real_t> &zvals) const;
  /**
   * @brief
   */
  void apply_initial_condition_xy();

  /**
   * @brief
   */
  void apply_initial_condition_xy(EulerBeamMesh &bmesh);

  /**
   * @brief
   */
  void initialize_quadrature_cache();

  /**
   * @brief
   */
  void compute_slopes_collocation();

  /**
   * @brief
   */
  void apply_initial_condition_pq();

  /**
   * @brief
   */
  void apply_boundary_condition_lambda();

  /**
   * @brief
   */
  virtual void update_pq();

  /**
   * @brief
   */
  void apply_boundary_condition_A();

  /**
   * @brief
   */
  void assemble_A();

  /**
   * @brief
   */
  void decompose_A();

  /**
   * @brief
   */
  void clear_rhs();

  /**
   * @brief
   */
  void assemble_constraint_rhs();

  /**
   * @brief
   */
  void add_uniform_load_rhs(std::array<real_t, 3> load);

  /**
   * @brief
   */
  void add_nodal_load_rhs(const std::vector<std::array<real_t, 3>> &load);

  /**
   * @brief
   */
  void add_point_boundary_loads();

  /**
   * @brief
   */
  void assemble_f(std::array<real_t, 3> load);

  /**
   * @brief
   */
  void apply_boundary_condition_f();

  /**
   * @brief
   */
  virtual void update_xy(std::array<real_t, 3> load);

  /**
   * @brief
   */
  virtual void update_multipliers();

  /**
   * @brief
   */
  bool is_converged(const VectorXd &x_old_iter, const VectorXd &x_new_iter,
                    const VectorXd &y_old_iter, const VectorXd &y_new_iter,
                    const VectorXd &z_old_iter,
                    const VectorXd &z_new_iter) const;

  /**
   * @brief
   */
  void assemble_mass_matrix();

  /**
   * @brief
   */
  void prepare_system_newmark(real_t dt);

  /**
   * @brief
   */
  void assemble_system_newmark_rhs(std::array<real_t, 3> load, real_t dt);

  /**
   * @brief
   */
  void
  assemble_system_newmark_rhs(const std::vector<std::array<real_t, 3>> &load,
                              real_t dt);

  /**
   * @brief
   */
  void assemble_system_newmark_rhs_averaged(
      const std::vector<std::array<real_t, 3>> &averaged_load, real_t dt);

  /**
   * @brief
   */
  void add_averaged_uniform_load_rhs(std::array<real_t, 3> load);

  /**
   * @brief
   */
  void
  add_averaged_nodal_load_rhs(const std::vector<std::array<real_t, 3>> &load);

  /**
   * @brief
   */
  void apply_midpoint_bending_rhs();

  /**
   * @brief
   */
  void update_average_acceleration_state_component(const VectorXd &u_old,
                                                   const VectorXd &u_new,
                                                   VectorXd &v_hist,
                                                   VectorXd &a_hist, real_t dt);
  /**
   * @brief
   */
  void apply_dynamic_state_boundary_conditions();

  /**
   * @brief
   */
  real_t compute_relative_state_update(const VectorXd &x_old_iter,
                                       const VectorXd &x_new_iter,
                                       const VectorXd &y_old_iter,
                                       const VectorXd &y_new_iter,
                                       const VectorXd &z_old_iter,
                                       const VectorXd &z_new_iter) const;
  /**
   * @brief
   */
  real_t compute_max_pq_error() const;

  /**
   * @brief
   */
  real_t compute_max_state_update(const VectorXd &x_old_iter,
                                  const VectorXd &x_new_iter,
                                  const VectorXd &y_old_iter,
                                  const VectorXd &y_new_iter,
                                  const VectorXd &z_old_iter,
                                  const VectorXd &z_new_iter) const;
};

} // namespace Models
} // namespace ELFF
