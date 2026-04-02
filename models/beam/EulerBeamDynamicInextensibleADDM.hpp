#pragma once

#include "elff/models/beam/EulerBeamStaticInextensibleADDM.hpp"

using namespace Eigen;

namespace ELFF {
namespace Models {

/**
 * @brief Dynamic inextensible Euler beam solved with the dense ADDM
 * formulation and the quasi-static time-discrete ALG2 formulation from
 * Basting et al.
 *
 * The implementation keeps the ADDM / ALG2 splitting from the static model
 * and applies it to the time-discrete problem at each step.  The public
 * `solve_newmark()` path matches the paper's average-acceleration Newmark
 * variant, where the inertia term acts on `x^{n+1}` while the bending and
 * load terms are midpoint-averaged.
 */
class EulerBeamDynamicInextensibleADDM : public EulerBeamStaticInextensibleADDM
{
public:
  EulerBeamDynamicInextensibleADDM(real_t length,
                                   real_t EI,
                                   real_t mu,
                                   size_t nodes,
                                   EulerBeam::EulerBeamBCs bcs,
                                   real_t r_penalty);

  void solve(real_t dt, std::array<real_t, 3> load) override;

  void solve(real_t dt, std::vector<std::array<real_t, 3>> load) override;

  void solve_newmark(real_t dt,
                     std::array<real_t, 3> load,
                     real_t beta,
                     real_t gamma);

  void solve_newmark(real_t dt,
                     std::vector<std::array<real_t, 3>> load,
                     real_t beta,
                     real_t gamma);

  void apply_initial_condition() override;

  void apply_initial_condition(EulerBeamMesh& bmesh) override;

protected:
  VectorXd x_prev;
  VectorXd y_prev;
  VectorXd z_prev;

  VectorXd vx_prev;
  VectorXd vy_prev;
  VectorXd vz_prev;

  VectorXd ax_prev;
  VectorXd ay_prev;
  VectorXd az_prev;

  MatrixXd mass;
  std::array<real_t, 3> load_prev;
  std::vector<std::array<real_t, 3>> nodal_load_prev;
  bool have_prev_uniform_load;
  bool have_prev_nodal_load;

  void assemble_mass_matrix();

  void assemble_system_newmark(std::array<real_t, 3> load, real_t dt);

  void assemble_system_newmark(const std::vector<std::array<real_t, 3>>& load,
                               real_t dt);

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

  void update_mesh();
};

} // namespace Models
} // namespace ELFF
