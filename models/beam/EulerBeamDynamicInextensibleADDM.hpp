#pragma once

#include "elff/models/beam/EulerBeamStaticInextensibleADDM.hpp"

using namespace Eigen;

namespace ELFF {
namespace Models {

/**
 * @brief Dynamic inextensible Euler beam solved with the dense ADDM
 * formulation and Newmark time integration.
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
  MatrixXd A_static_unconstrained;
  std::array<real_t, 3> load_prev;

  void assemble_mass_matrix();

  void assemble_system_newmark(std::array<real_t, 3> load,
                               real_t dt,
                               real_t beta,
                               real_t gamma);

  void assemble_system_newmark(
    const std::vector<std::array<real_t, 3>>& load,
    real_t dt,
    real_t beta,
    real_t gamma);

  void assemble_f_nodal(const std::vector<std::array<real_t, 3>>& load);

  void initialize_newmark_acceleration(std::array<real_t, 3> load);

  void initialize_newmark_acceleration(
    const std::vector<std::array<real_t, 3>>& load);

  void apply_dynamic_state_boundary_conditions();

  real_t compute_relative_xy_update(const VectorXd& x_old_iter,
                                    const VectorXd& x_new_iter,
                                    const VectorXd& y_old_iter,
                                    const VectorXd& y_new_iter) const;

  real_t compute_max_pq_error() const;

  real_t compute_max_xy_update(const VectorXd& x_old_iter,
                               const VectorXd& x_new_iter,
                               const VectorXd& y_old_iter,
                               const VectorXd& y_new_iter) const;

  void update_newmark_state_component(const VectorXd& u_old,
                                      const VectorXd& u_new,
                                      VectorXd& v_hist,
                                      VectorXd& a_hist,
                                      real_t dt,
                                      real_t beta,
                                      real_t gamma);

  void update_mesh();
};

} // namespace Models
} // namespace ELFF
