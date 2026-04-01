#pragma once

#include "elff/models/beam/EulerBeamStaticInextensibleAugKKT.hpp"

using namespace Eigen;

namespace ELFF {
namespace Models {

class EulerBeamDynamicInextensibleAugKKT
  : public EulerBeamStaticInextensibleAugKKT
{
public:
  EulerBeamDynamicInextensibleAugKKT(real_t length,
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
  VectorXd v_prev;
  VectorXd a_prev;
  VectorXd u_prev;
  SparseMatrix<real_t> mass;
  std::array<real_t, 3> load_prev;

  void assemble_mass_matrix();

  void assemble_system_newmark(real_t dt,
                               std::array<real_t, 3> load,
                               real_t beta,
                               real_t gamma);

  void assemble_system_newmark(real_t dt,
                               std::vector<std::array<real_t, 3>> load,
                               real_t beta,
                               real_t gamma);

  void assemble_system(std::vector<std::array<real_t, 3>> load);

  void add_newmark_inertial_terms(real_t dt, real_t beta);

  void apply_dynamic_state_boundary_conditions();

  void update_mesh();
};

} // namespace Models
} // namespace ELFF
