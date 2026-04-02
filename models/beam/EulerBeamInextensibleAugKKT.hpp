#pragma once

#include "elff/fem/Shapes.hpp"
#include "elff/models/beam/EulerBeam.hpp"

#include <Eigen/Dense>
#include <Eigen/SparseLU>
#include <unsupported/Eigen/AutoDiff>

using namespace Eigen;

namespace ELFF {
namespace Models {

/**
 * @brief Inextensible Euler beam solved with a sparse augmented-KKT
 * formulation.
 *
 * The class supports both static and dynamic solves. Static instances are
 * constructed without a mass density, while dynamic instances use Newmark time
 * integration.
 */
class EulerBeamInextensibleAugKKT : public EulerBeam
{
public:
  EulerBeamInextensibleAugKKT(real_t length,
                              real_t EI,
                              size_t nodes,
                              EulerBeam::EulerBeamBCs bcs,
                              real_t r_penalty);

  EulerBeamInextensibleAugKKT(real_t length,
                              real_t EI,
                              real_t mu,
                              size_t nodes,
                              EulerBeam::EulerBeamBCs bcs,
                              real_t r_penalty);

  ~EulerBeamInextensibleAugKKT();

  void solve() override;

  void solve(std::array<real_t, 3> load) override;

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

  virtual void apply_initial_condition();

  virtual void apply_initial_condition(EulerBeamMesh& bmesh) override;

protected:
  size_t dimension;
  size_t elements, nodes;
  real_t ds;
  size_t ndof_x, ndof_y, ndof_z, ndof_l;
  size_t offset_x, offset_y, offset_z, offset_l;
  size_t ndof;
  real_t r_penalty;
  size_t max_iter;
  real_t tol;

  VectorXd residual;
  SparseMatrix<real_t> jacobian;
  VectorXd u;

  VectorXd v_prev;
  VectorXd a_prev;
  VectorXd u_prev;
  SparseMatrix<real_t> mass;
  std::array<real_t, 3> load_prev;

  void assemble_residual(std::array<real_t, 3> load);

  void assemble_system(std::array<real_t, 3> load);

  void assemble_system(std::vector<std::array<real_t, 3>> load);

  void assemble_system_newmark(real_t dt,
                               std::array<real_t, 3> load,
                               real_t beta,
                               real_t gamma);

  void assemble_system_newmark(real_t dt,
                               std::vector<std::array<real_t, 3>> load,
                               real_t beta,
                               real_t gamma);

  void assemble_mass_matrix();

  void add_newmark_inertial_terms(real_t dt, real_t beta);

  void apply_boundary_conditions();

  void apply_dynamic_state_boundary_conditions();

  void update_mesh();

  std::array<size_t, 14> get_element_dof_indices(size_t e) const;

  Matrix<real_t, 14, 1> get_element_state(
    const std::array<size_t, 14>& idx) const;

  template<typename T>
  Matrix<T, 14, 1> assemble_element_residual_template(
    const Matrix<T, 14, 1>& u_elem,
    const std::array<real_t, 3> load) const
  {
    Matrix<T, 14, 1> residual = Matrix<T, 14, 1>::Zero();

    const Matrix<T, 4, 1> ux = u_elem.template segment<4>(0);
    const Matrix<T, 4, 1> uy = u_elem.template segment<4>(4);
    const Matrix<T, 4, 1> uz = u_elem.template segment<4>(8);
    const Matrix<T, 2, 1> ul = u_elem.template segment<2>(12);

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

      const T l = M[0] * ul(0) + M[1] * ul(1);
      const T S = xp * xp + yp * yp + zp * zp - 1.0;

      for (size_t a = 0; a < 4; ++a) {
        residual(a) += EI * xpp * ddH[a] * w * ds;
        residual(4 + a) += EI * ypp * ddH[a] * w * ds;
        residual(8 + a) += EI * zpp * ddH[a] * w * ds;

        residual(a) -= load[0] * H[a] * w * ds;
        residual(4 + a) -= load[1] * H[a] * w * ds;
        residual(8 + a) -= load[2] * H[a] * w * ds;

        const T coeff = 2 * (l + r_penalty * S) * dH[a] * w * ds;
        residual(a) += xp * coeff;
        residual(4 + a) += yp * coeff;
        residual(8 + a) += zp * coeff;
      }

      for (size_t a = 0; a < 2; ++a) {
        residual(12 + a) += S * M[a] * w * ds;
      }
    }

    return residual;
  }

  template<typename T>
  Matrix<T, 14, 1> assemble_element_residual_template(
    const Matrix<T, 14, 1>& u_elem,
    const std::array<std::array<real_t, 3>, 2>& load_elem) const
  {
    Matrix<T, 14, 1> residual = Matrix<T, 14, 1>::Zero();

    const Matrix<T, 4, 1> ux = u_elem.template segment<4>(0);
    const Matrix<T, 4, 1> uy = u_elem.template segment<4>(4);
    const Matrix<T, 4, 1> uz = u_elem.template segment<4>(8);
    const Matrix<T, 2, 1> ul = u_elem.template segment<2>(12);

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

      const T l = M[0] * ul(0) + M[1] * ul(1);
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

        const T coeff = 2 * (l + r_penalty * S) * dH[a] * w * ds;
        residual(a) += xp * coeff;
        residual(4 + a) += yp * coeff;
        residual(8 + a) += zp * coeff;
      }

      for (size_t a = 0; a < 2; ++a) {
        residual(12 + a) += S * M[a] * w * ds;
      }
    }

    return residual;
  }
};

} // namespace Models
} // namespace ELFF
