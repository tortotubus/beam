#pragma once

#include "elff/fem/Shapes.hpp"
#include "elff/models/beam/EulerBeam.hpp"

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
 * @brief Inextensible Euler beam solved with a sparse finite-element penalty
 * formulation.
 *
 * The same finite-element residual is used for both static equilibrium solves
 * and dynamic Newmark time stepping. Inextensibility is enforced through a
 * quadratic penalty on the strain defect
 *
 *   S(u) = |r'(s)|^2 - 1.
 */
class EulerBeamInextensiblePenalty : public EulerBeam
{
public:
  EulerBeamInextensiblePenalty(real_t                  length,
                               real_t                  EI,
                               size_t                  nodes,
                               EulerBeam::EulerBeamBCs bcs,
                               real_t                  r_penalty);

  EulerBeamInextensiblePenalty(real_t                  length,
                               real_t                  EI,
                               real_t                  mu,
                               size_t                  nodes,
                               EulerBeam::EulerBeamBCs bcs,
                               real_t                  r_penalty);

  ~EulerBeamInextensiblePenalty() = default;

  void solve() override;
  void solve(std::array<real_t, 3> load) override;
  void solve(std::vector<std::array<real_t, 3>> load) override;
  void solve(real_t dt, std::array<real_t, 3> load) override;
  void solve(real_t dt, std::vector<std::array<real_t, 3>> load) override;

  void apply_initial_condition();
  void apply_initial_condition(EulerBeamMesh& bmesh) override;

protected:
  size_t dimension;
  size_t elements;
  size_t nodes;
  real_t ds;

  size_t ndof_x;
  size_t ndof_y;
  size_t ndof_z;
  size_t offset_x;
  size_t offset_y;
  size_t offset_z;
  size_t ndof;

  real_t r_penalty;
  size_t max_iter_nonlinear;
  real_t tol_primal;

  VectorXd             residual;
  SparseMatrix<real_t> jacobian;
  VectorXd             u;
  VectorXd             u_prev;
  VectorXd             v_prev;
  VectorXd             a_prev;

  std::array<size_t, 12> get_element_dof_indices(size_t e) const;
  Matrix<real_t, 12, 1> get_element_state(
    const std::array<size_t, 12>& idx) const;

  template<typename T>
  Matrix<T, 12, 1> assemble_element_residual_template(
    const Matrix<T, 12, 1>&       u_elem,
    const std::array<real_t, 3>& load) const
  {
    Matrix<T, 12, 1> residual = Matrix<T, 12, 1>::Zero();

    const Matrix<T, 4, 1> ux = u_elem.template segment<4>(0);
    const Matrix<T, 4, 1> uy = u_elem.template segment<4>(4);
    const Matrix<T, 4, 1> uz = u_elem.template segment<4>(8);

    const real_t xi_q[] = { 0.1127016654, 0.5, 0.8872983346 };
    const real_t w_q[] = { 0.2777777778, 0.4444444444, 0.2777777778 };

    for (size_t qi = 0; qi < 3; ++qi) {
      const real_t xi = xi_q[qi];
      const real_t w = w_q[qi];

      const auto H = ELFF::FEM::CubicHermite<real_t>::values(xi, ds);
      const auto dH = ELFF::FEM::CubicHermite<real_t>::derivs(xi, ds);
      const auto ddH = ELFF::FEM::CubicHermite<real_t>::second_derivs(xi, ds);

      T xp = 0, xpp = 0;
      T yp = 0, ypp = 0;
      T zp = 0, zpp = 0;

      for (size_t i = 0; i < 4; ++i) {
        xp += dH[i] * ux(i);
        yp += dH[i] * uy(i);
        zp += dH[i] * uz(i);
        xpp += ddH[i] * ux(i);
        ypp += ddH[i] * uy(i);
        zpp += ddH[i] * uz(i);
      }

      const T S = xp * xp + yp * yp + zp * zp - 1.0;

      for (size_t a = 0; a < 4; ++a) {
        residual(a) += EI * xpp * ddH[a] * w * ds;
        residual(4 + a) += EI * ypp * ddH[a] * w * ds;
        residual(8 + a) += EI * zpp * ddH[a] * w * ds;

        residual(a) -= load[0] * H[a] * w * ds;
        residual(4 + a) -= load[1] * H[a] * w * ds;
        residual(8 + a) -= load[2] * H[a] * w * ds;

        const T coeff = 2.0 * r_penalty * S * dH[a] * w * ds;
        residual(a) += xp * coeff;
        residual(4 + a) += yp * coeff;
        residual(8 + a) += zp * coeff;
      }
    }

    return residual;
  }

  template<typename T>
  Matrix<T, 12, 1> assemble_element_residual_template(
    const Matrix<T, 12, 1>&                     u_elem,
    const std::array<std::array<real_t, 3>, 2>& load_elem) const
  {
    Matrix<T, 12, 1> residual = Matrix<T, 12, 1>::Zero();

    const Matrix<T, 4, 1> ux = u_elem.template segment<4>(0);
    const Matrix<T, 4, 1> uy = u_elem.template segment<4>(4);
    const Matrix<T, 4, 1> uz = u_elem.template segment<4>(8);

    const real_t xi_q[] = { 0.1127016654, 0.5, 0.8872983346 };
    const real_t w_q[] = { 0.2777777778, 0.4444444444, 0.2777777778 };

    for (size_t qi = 0; qi < 3; ++qi) {
      const real_t xi = xi_q[qi];
      const real_t w = w_q[qi];

      const auto H = ELFF::FEM::CubicHermite<real_t>::values(xi, ds);
      const auto dH = ELFF::FEM::CubicHermite<real_t>::derivs(xi, ds);
      const auto ddH = ELFF::FEM::CubicHermite<real_t>::second_derivs(xi, ds);
      const auto M = ELFF::FEM::LinearShape<real_t>::values(xi);

      T xp = 0, xpp = 0;
      T yp = 0, ypp = 0;
      T zp = 0, zpp = 0;

      for (size_t i = 0; i < 4; ++i) {
        xp += dH[i] * ux(i);
        yp += dH[i] * uy(i);
        zp += dH[i] * uz(i);
        xpp += ddH[i] * ux(i);
        ypp += ddH[i] * uy(i);
        zpp += ddH[i] * uz(i);
      }

      const real_t fx = M[0] * load_elem[0][0] + M[1] * load_elem[1][0];
      const real_t fy = M[0] * load_elem[0][1] + M[1] * load_elem[1][1];
      const real_t fz = M[0] * load_elem[0][2] + M[1] * load_elem[1][2];
      const T S = xp * xp + yp * yp + zp * zp - 1.0;

      for (size_t a = 0; a < 4; ++a) {
        residual(a) += EI * xpp * ddH[a] * w * ds;
        residual(4 + a) += EI * ypp * ddH[a] * w * ds;
        residual(8 + a) += EI * zpp * ddH[a] * w * ds;

        residual(a) -= fx * H[a] * w * ds;
        residual(4 + a) -= fy * H[a] * w * ds;
        residual(8 + a) -= fz * H[a] * w * ds;

        const T coeff = 2.0 * r_penalty * S * dH[a] * w * ds;
        residual(a) += xp * coeff;
        residual(4 + a) += yp * coeff;
        residual(8 + a) += zp * coeff;
      }
    }

    return residual;
  }

  void assemble_system(std::array<real_t, 3> load);
  void assemble_system(std::vector<std::array<real_t, 3>> load);

  void apply_boundary_conditions();
  void update_mesh();

  void solve_newmark(real_t dt,
                     std::array<real_t, 3> load,
                     real_t beta,
                     real_t gamma);
  void solve_newmark(real_t dt,
                     std::vector<std::array<real_t, 3>> load,
                     real_t beta,
                     real_t gamma);

  void initialize_acceleration(std::array<real_t, 3> load);
  void initialize_acceleration(std::vector<std::array<real_t, 3>> load);
  void add_newmark_inertial_terms(real_t dt, real_t beta);
  void update_newmark_state_from_displacement(real_t dt,
                                              real_t beta,
                                              real_t gamma,
                                              const VectorXd& v_old,
                                              const VectorXd& a_old);
  void apply_dynamic_state_boundary_conditions();
};

} // namespace Models
} // namespace ELFF
