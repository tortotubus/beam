#pragma once

#include "elff/fem/Shapes.hpp"
#include "elff/models/beam/EulerBeam.hpp"

#include <Eigen/Dense>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseLU>

using namespace Eigen;

namespace ELFF {
namespace Models {

/**
 * @brief Dynamic inextensible Euler beam solved as an index-1 DAE with a
 * sparse saddle-point system and Newmark time integration.
 *
 * This class is a standalone dynamic beam implementation. Static solve
 * overloads are intentionally left unimplemented.
 */
class EulerBeamInextensibleIndex1 : public EulerBeam
{
public:
  EulerBeamInextensibleIndex1(real_t                  length,
                                     real_t                  EI,
                                     real_t                  mu,
                                     size_t                  n_nodes,
                                     EulerBeam::EulerBeamBCs bcs);

  void solve() override;
  void solve(std::array<real_t, 3> load) override;
  void solve(std::vector<std::array<real_t, 3>> load) override;

  void apply_initial_condition();
  void apply_initial_condition(EulerBeamMesh& bmesh) override;

  void solve(real_t dt, std::array<real_t, 3> load) override;
  void solve(real_t dt, std::vector<std::array<real_t, 3>> load) override;

protected:
  size_t elements, nodes;
  real_t ds;
  size_t ndof_x, ndof_y, ndof_z, ndof_l;
  size_t offset_x, offset_y, offset_z, offset_l;
  size_t ndof;

  VectorXd u;
  VectorXd lambda;
  VectorXd u_prev;
  VectorXd v_prev;
  VectorXd a_prev;

  SparseMatrix<real_t> K_elastic;
  SparseMatrix<real_t> saddle_mat;
  VectorXd saddle_rhs;

  real_t tol_position_drift;
  real_t tol_velocity_drift;
  size_t max_projection_iter;

  void assemble_elastic_stiffness();
  void assemble_B_and_g(SparseMatrix<real_t>& B, VectorXd& g) const;
  VectorXd assemble_constraint_residual_vector() const;
  VectorXd assemble_f_ext(std::array<real_t, 3> load) const;
  VectorXd assemble_f_ext(const std::vector<std::array<real_t, 3>>& load) const;
  void assemble_saddle_system(real_t                      dt,
                              real_t                      beta,
                              const VectorXd&             f_ext,
                              const SparseMatrix<real_t>& B,
                              const VectorXd&             g);
  void apply_saddle_boundary_conditions();
  void update_newmark_state(real_t dt, real_t beta, real_t gamma);
  void apply_dynamic_state_boundary_conditions();
  void project_position_onto_constraint();
  void project_velocity_onto_constraint(const SparseMatrix<real_t>& B);
  void update_mesh();
};

} // namespace Models
} // namespace ELFF
