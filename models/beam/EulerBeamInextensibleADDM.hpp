#pragma once

#include "elff/fem/Shapes.hpp"
#include "elff/models/beam/EulerBeam.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <Eigen/IterativeLinearSolvers>
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
class EulerBeamInextensibleADDM : public EulerBeam
{
public:
  EulerBeamInextensibleADDM(real_t length,
                            real_t EI,
                            size_t nodes,
                            EulerBeamBCs bcs,
                            real_t r_penalty);

  EulerBeamInextensibleADDM(real_t length,
                            real_t EI,
                            real_t mu,
                            size_t nodes,
                            EulerBeamBCs bcs,
                            real_t r_penalty);

  ~EulerBeamInextensibleADDM();

  virtual void solve() override;
  virtual void solve(std::array<real_t, 3> load) override;
  virtual void solve(real_t dt, std::array<real_t, 3> load) override;
  virtual void solve(real_t dt, std::vector<std::array<real_t, 3>> load) override;

  void solve_newmark(real_t dt,
                     std::array<real_t, 3> load,
                     real_t beta,
                     real_t gamma);

  void solve_newmark(real_t dt,
                     std::vector<std::array<real_t, 3>> load,
                     real_t beta,
                     real_t gamma);

  virtual void apply_initial_condition();
  virtual void apply_initial_condition(EulerBeamMesh& bmesh) override;

  const VectorXd& get_lambda_x() const;
  const VectorXd& get_lambda_y() const;
  const VectorXd& get_lambda_z() const;
  const MatrixXd& get_A() const;

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

  VectorXd lambda_x, lambda_y, lambda_z;
  VectorXd p, q, r;
  VectorXd xp, yp, zp;

  size_t max_outer;
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
  VectorXd mass_diag;
  std::array<real_t, 3> load_prev;
  std::vector<std::array<real_t, 3>> nodal_load_prev;
  bool have_prev_uniform_load;
  bool have_prev_nodal_load;

  void update_mesh();
  void collect_boundary_dofs(std::vector<size_t>& idx,
                              std::vector<real_t>& xvals,
                              std::vector<real_t>& yvals,
                              std::vector<real_t>& zvals) const;
  void apply_initial_condition_xy();
  void apply_initial_condition_xy(EulerBeamMesh& bmesh);
  void initialize_quadrature_cache();
  void compute_slopes_collocation();
  void apply_initial_condition_pq();
  void apply_boundary_condition_lambda();
  virtual void update_pq();
  void apply_boundary_condition_A();
  void assemble_A();
  void decompose_A();
  void clear_rhs();
  void assemble_constraint_rhs();
  void add_uniform_load_rhs(std::array<real_t, 3> load);
  void add_nodal_load_rhs(const std::vector<std::array<real_t, 3>>& load);
  void add_point_boundary_loads();
  void assemble_f(std::array<real_t, 3> load);
  void apply_boundary_condition_f();
  virtual void update_xy(std::array<real_t, 3> load);
  virtual void update_multipliers();
  bool is_converged(const VectorXd& x_old_iter,
                    const VectorXd& x_new_iter,
                    const VectorXd& y_old_iter,
                    const VectorXd& y_new_iter,
                    const VectorXd& z_old_iter,
                    const VectorXd& z_new_iter) const;

  void assemble_mass_matrix();
  void prepare_system_newmark(real_t dt);
  void assemble_system_newmark_rhs(std::array<real_t, 3> load, real_t dt);
  void assemble_system_newmark_rhs(
    const std::vector<std::array<real_t, 3>>& load, real_t dt);
  void add_averaged_uniform_load_rhs(std::array<real_t, 3> load);
  void add_averaged_nodal_load_rhs(
    const std::vector<std::array<real_t, 3>>& load);
  void apply_midpoint_bending_rhs();
  void update_average_acceleration_state_component(const VectorXd& u_old,
                                                   const VectorXd& u_new,
                                                   VectorXd& v_hist,
                                                   VectorXd& a_hist,
                                                   real_t dt);
  void apply_dynamic_state_boundary_conditions();
  real_t compute_relative_state_update(const VectorXd& x_old_iter,
                                       const VectorXd& x_new_iter,
                                       const VectorXd& y_old_iter,
                                       const VectorXd& y_new_iter,
                                       const VectorXd& z_old_iter,
                                       const VectorXd& z_new_iter) const;
  real_t compute_max_pq_error() const;
  real_t compute_max_state_update(const VectorXd& x_old_iter,
                                  const VectorXd& x_new_iter,
                                  const VectorXd& y_old_iter,
                                  const VectorXd& y_new_iter,
                                  const VectorXd& z_old_iter,
                                  const VectorXd& z_new_iter) const;
};

} // namespace Models
} // namespace ELFF
