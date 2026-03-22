#pragma once

// #include <beam/LinAlg/Matrix.hpp>
// #include <beam/LinAlg/Vector.hpp>

#include "EulerBeam.hpp"
#include "Shapes.hpp"

#include <cstdio> // for popen, pclose, fprintf
#include <iostream>
#include <stdlib.h>
#include <utility>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCholesky> // for SimplicialLLT
#include <unsupported/Eigen/AutoDiff>

using namespace Eigen;

namespace ELFF {
namespace Models {

class EulerBeamStaticInextensibleMoMSparse : public EulerBeam
{
public:
  EulerBeamStaticInextensibleMoMSparse(real_t length,
                                       real_t EI,
                                       size_t nodes,
                                       EulerBeam::EulerBeamBCs bcs,
                                       real_t r_penalty)
    : EulerBeam(length, EI, nodes, bcs)
    , elements(nodes - 1)
    , nodes(nodes)
    , ds(mesh.get_ds())
    , ndof_x(2 * nodes)
    , ndof_y(2 * nodes)
    , ndof_z(2 * nodes)
    , ndof_l(1 * nodes)
    , ndof(ndof_x + ndof_y + ndof_z)
    , offset_x(0)
    , offset_y(ndof_x)
    , offset_z(ndof_x + ndof_y)
    // , offset_l(ndof_x + ndof_y + ndof_z)
    , jacobian(SparseMatrix<real_t>(ndof, ndof))
    , residual(VectorXd::Zero(ndof))
    , lambda(VectorXd::Zero(ndof_l))
    , u(VectorXd::Zero(ndof))
    , r_penalty(r_penalty)
    , max_iter_inner(1000)
    , max_iter_outer(1000)
    , tol_inner(1e-5)
    , tol_outer(1e-5)
  {
    apply_initial_condition(mesh);
  };

  ~EulerBeamStaticInextensibleMoMSparse() {};

  void solve() override { solve({ 0., 0., 0. }); }

  void solve(std::array<real_t, 3> load) override
  {
    real_t S_norm = 0;

    // ConjugateGradient<
    //   SparseMatrix<real_t>,          // or SparseMatrix<double>
    //   Lower | Upper,          // tell it K is symmetric
    //   DiagonalPreconditioner<real_t> // Jacobi preconditioner
    // >
    //   solver;

    // SimplicialLLT<SparseMatrix<real_t>> solver;

    ConjugateGradient<SparseMatrix<real_t>,
                      Lower | Upper,
                      IncompleteCholesky<real_t>>
      solver;

    for (size_t iter_outer = 0; iter_outer < max_iter_outer; iter_outer++) {
      assemble_system(load);
      apply_boundary_conditions();

      real_t res_norm = residual.norm();

      if (res_norm < tol_outer) {
        break;
      } else if (iter_outer == max_iter_outer - 1) {
        ELFF_ABORT(
          "EulerBeamStaticInexntensibleMoM::solve() did not converge.\n");
      }

      solver.setTolerance(tol_inner);
      // solver.setMaxIterations(100);
      solver.compute(jacobian);

      if (solver.info() != Success) {
        ELFF_ABORT("EulerBeamStaticInextensibleMoMSparse::solve(): "
                   "Preconditioner failed.\n");
      }

      VectorXd delta_u = solver.solve(-residual);
      u += delta_u;

      S_norm = update_lambda();
    }

    update_mesh();
  }

  void solve(std::vector<std::array<real_t, 3>> load) override
  {
    ELFF_ASSERT(load.size() == nodes,
                "Size of load vector must equal number of nodes.");

    real_t S_norm = 0;

    ConjugateGradient<SparseMatrix<real_t>,
                      Lower | Upper,
                      IncompleteCholesky<real_t>>
      solver;

    for (size_t iter_outer = 0; iter_outer < max_iter_outer; iter_outer++) {
      assemble_system(load);
      apply_boundary_conditions();

      real_t res_norm = residual.norm();

      if (res_norm < tol_outer) {
        break;
      } else if (iter_outer == max_iter_outer - 1) {
        ELFF_ABORT(
          "EulerBeamStaticInexntensibleMoM::solve() did not converge.\n");
      }

      solver.setTolerance(tol_inner);
      solver.compute(jacobian);

      if (solver.info() != Success) {
        ELFF_ABORT("EulerBeamStaticInextensibleMoMSparse::solve(): "
                   "Preconditioner failed.\n");
      }

      VectorXd delta_u = solver.solve(-residual);
      u += delta_u;

      S_norm = update_lambda();
    }

    (void) S_norm;

    update_mesh();
  }

  virtual void apply_initial_condition(EulerBeamMesh& bmesh) override
  {
    ELFF_ASSERT(
      nodes == bmesh.get_nodes(),
      "Provided mesh must have same node count as the previous mesh.\n");

    auto centerline = bmesh.get_centerline();
    auto slopes = bmesh.get_slope();

    for (size_t ni = 0; ni < nodes; ni++) {
      u(offset_x + 2 * ni + 0) = centerline[ni][0];
      u(offset_x + 2 * ni + 1) = slopes[ni][0];
      u(offset_y + 2 * ni + 0) = centerline[ni][1];
      u(offset_y + 2 * ni + 1) = slopes[ni][1];
      u(offset_z + 2 * ni + 0) = centerline[ni][2];
      u(offset_z + 2 * ni + 1) = slopes[ni][2];
    }

    update_mesh();
  }

  virtual void apply_initial_condition()
  {
    for (size_t i = 0; i < nodes; i++) {
      u(offset_x + 2 * i + 0) = ds * i;
      u(offset_x + 2 * i + 1) = 1.;
      u(offset_y + 2 * i + 0) = 0.;
      u(offset_y + 2 * i + 1) = 0.;
      u(offset_z + 2 * i + 0) = 0.;
      u(offset_z + 2 * i + 1) = 0.;
    }
  }

protected:
  size_t dimension;
  size_t elements, nodes;
  real_t ds;
  size_t ndof_x, ndof_y, ndof_z, ndof_l;
  size_t offset_x, offset_y, offset_z, offset_l;
  size_t ndof;
  real_t r_penalty;
  size_t max_iter_inner, max_iter_outer;
  real_t tol_inner, tol_outer;

  VectorXd residual, lambda;
  SparseMatrix<real_t> jacobian;
  VectorXd u;

  EulerBeamStaticInextensibleMoMSparse(real_t length,
                                       real_t EI,
                                       real_t mu,
                                       size_t nodes,
                                       EulerBeam::EulerBeamBCs bcs,
                                       real_t r_penalty)
    : EulerBeam(length, EI, mu, nodes, bcs)
    , elements(nodes - 1)
    , nodes(nodes)
    , ds(mesh.get_ds())
    , ndof_x(2 * nodes)
    , ndof_y(2 * nodes)
    , ndof_z(2 * nodes)
    , ndof_l(1 * nodes)
    , ndof(ndof_x + ndof_y + ndof_z)
    , offset_x(0)
    , offset_y(ndof_x)
    , offset_z(ndof_x + ndof_y)
    , jacobian(SparseMatrix<real_t>(ndof, ndof))
    , residual(VectorXd::Zero(ndof))
    , lambda(VectorXd::Zero(ndof_l))
    , u(VectorXd::Zero(ndof))
    , r_penalty(r_penalty)
    , max_iter_inner(1000)
    , max_iter_outer(1000)
    , tol_inner(1e-6)
    , tol_outer(1e-6)
  {
    apply_initial_condition(mesh);
  };

  void assemble_residual(std::array<real_t, 3> load)
  {
    residual = assemble_residual_template<real_t>(u, load);
  }

  void assemble_residual(std::vector<std::array<real_t, 3>> load)
  {
    residual = assemble_residual_template<real_t>(u, load);
  }

  std::array<size_t, 12> get_element_dof_indices(size_t e) const
  {
    const size_t n0 = e;
    const size_t n1 = e + 1;

    return { offset_x + 2 * n0 + 0,
             offset_x + 2 * n0 + 1,
             offset_x + 2 * n1 + 0,
             offset_x + 2 * n1 + 1,
             offset_y + 2 * n0 + 0,
             offset_y + 2 * n0 + 1,
             offset_y + 2 * n1 + 0,
             offset_y + 2 * n1 + 1,
             offset_z + 2 * n0 + 0,
             offset_z + 2 * n0 + 1,
             offset_z + 2 * n1 + 0,
             offset_z + 2 * n1 + 1 };
  }

  Matrix<real_t, 12, 1> get_element_state(
    const std::array<size_t, 12>& idx) const
  {
    Matrix<real_t, 12, 1> u_elem;
    for (int i = 0; i < 12; ++i) {
      u_elem(i) = u(idx[i]);
    }
    return u_elem;
  }

  std::array<real_t, 2> get_element_lambda(size_t e) const
  {
    return { lambda(e), lambda(e + 1) };
  }

  template<typename T>
  Matrix<T, 12, 1> assemble_element_residual_template(
    const Matrix<T, 12, 1>& u_elem,
    const std::array<real_t, 2>& lambda_elem,
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

      const auto H = CubicHermite<real_t>::values(xi, ds);
      const auto dH = CubicHermite<real_t>::derivs(xi, ds);
      const auto ddH = CubicHermite<real_t>::second_derivs(xi, ds);
      const auto M = LinearShape<real_t>::values(xi);

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

      const T l = M[0] * lambda_elem[0] + M[1] * lambda_elem[1];
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
    }

    return residual;
  }

  template<typename T>
  Matrix<T, 12, 1> assemble_element_residual_template(
    const Matrix<T, 12, 1>& u_elem,
    const std::array<real_t, 2>& lambda_elem,
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

      const auto H = CubicHermite<real_t>::values(xi, ds);
      const auto dH = CubicHermite<real_t>::derivs(xi, ds);
      const auto ddH = CubicHermite<real_t>::second_derivs(xi, ds);
      const auto M = LinearShape<real_t>::values(xi);

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

      const T l = M[0] * lambda_elem[0] + M[1] * lambda_elem[1];
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
    }

    return residual;
  }

  /**
   *
   */
  real_t update_lambda(real_t omega = 1.0)
  {
    real_t xi_q[] = { 0.1127016654, 0.5, 0.8872983346 };
    real_t w_q[] = { 0.2777777778, 0.4444444444, 0.2777777778 };

    Matrix<real_t, Dynamic, 1> lambda_n =
      Matrix<real_t, Dynamic, 1>::Zero(ndof_l);

    for (size_t e = 0; e < elements; ++e) {
      std::vector<size_t> elem_nodes = { e, e + 1 };
      std::vector<size_t> idx_x = { offset_x + 2 * elem_nodes[0],
                                    offset_x + 2 * elem_nodes[0] + 1,
                                    offset_x + 2 * elem_nodes[1],
                                    offset_x + 2 * elem_nodes[1] + 1 };
      std::vector<size_t> idx_y = { offset_y + 2 * elem_nodes[0],
                                    offset_y + 2 * elem_nodes[0] + 1,
                                    offset_y + 2 * elem_nodes[1],
                                    offset_y + 2 * elem_nodes[1] + 1 };
      std::vector<size_t> idx_z = { offset_z + 2 * elem_nodes[0],
                                    offset_z + 2 * elem_nodes[0] + 1,
                                    offset_z + 2 * elem_nodes[1],
                                    offset_z + 2 * elem_nodes[1] + 1 };
      std::vector<size_t> idx_l = { elem_nodes[0], elem_nodes[1] };

      std::array<real_t, 4> ux = {
        u[idx_x[0]], u[idx_x[1]], u[idx_x[2]], u[idx_x[3]]
      };
      std::array<real_t, 4> uy = {
        u[idx_y[0]], u[idx_y[1]], u[idx_y[2]], u[idx_y[3]]
      };
      std::array<real_t, 4> uz = {
        u[idx_z[0]], u[idx_z[1]], u[idx_z[2]], u[idx_z[3]]
      };
      std::array<real_t, 2> ul = { lambda[idx_l[0]], lambda[idx_l[1]] };

      std::vector<real_t> R_loc_l(2, 0);

      for (size_t qi = 0; qi < 3; ++qi) {
        real_t xi = xi_q[qi];
        real_t w = w_q[qi];

        auto H = CubicHermite<real_t>::values(xi, ds);
        auto dH = CubicHermite<real_t>::derivs(xi, ds);
        auto ddH = CubicHermite<real_t>::second_derivs(xi, ds);
        auto M = LinearShape<real_t>::values(xi);

        real_t x = 0, xp = 0, xpp = 0;
        real_t y = 0, yp = 0, ypp = 0;
        real_t z = 0, zp = 0, zpp = 0;

        for (size_t i = 0; i < 4; i++) {
          x += H[i] * ux[i];
          xp += dH[i] * ux[i];
          xpp += ddH[i] * ux[i];
          y += H[i] * uy[i];
          yp += dH[i] * uy[i];
          ypp += ddH[i] * uy[i];
          z += H[i] * uz[i];
          zp += dH[i] * uz[i];
          zpp += ddH[i] * uz[i];
        }

        real_t l = 0;
        for (size_t i = 0; i < 2; i++) {
          l += M[i] * ul[i];
        }

        real_t S = xp * xp + yp * yp + zp * zp - 1.0;

        for (size_t a = 0; a < 2; ++a) {
          // Variation w.r.t. lambda
          R_loc_l[a] += S * M[a] * w * ds;
        }
      }

      for (size_t i = 0; i < 2; ++i) {
        lambda_n[idx_l[i]] += R_loc_l[i];
      }
    }

    lambda = lambda_n;

    return lambda_n.norm();
  }

  void update_mesh()
  {
    size_t nodes = this->mesh.get_nodes();

    std::vector<std::array<real_t, 3>>& centerline = mesh.get_centerline();
    std::vector<std::array<real_t, 3>>& slope = mesh.get_slope();
    std::vector<real_t>& s = mesh.get_curvilinear_axis();

    for (size_t i = 0; i < nodes; ++i) {
      centerline[i][0] = u(offset_x + 2 * i);
      centerline[i][1] = u(offset_y + 2 * i);
      centerline[i][2] = u(offset_z + 2 * i);
      slope[i][0] = u(offset_x + 2 * i + 1);
      slope[i][1] = u(offset_y + 2 * i + 1);
      slope[i][2] = u(offset_z + 2 * i + 1);
    }
  }

  void assemble_system(std::array<real_t, 3> load)
  {
    using ADDeriv = Matrix<real_t, 12, 1>;
    using AD = AutoDiffScalar<ADDeriv>;
    using ADVec = Matrix<AD, 12, 1>;
    using Tpl = Triplet<real_t>;

    residual = VectorXd::Zero(ndof);
    std::vector<Tpl> triplets;
    triplets.reserve(elements * 12 * 12);

    for (size_t e = 0; e < elements; ++e) {
      const auto idx = get_element_dof_indices(e);
      const auto lambda_elem = get_element_lambda(e);
      const Matrix<real_t, 12, 1> u_elem = get_element_state(idx);

      ADVec u_ad;
      for (int a = 0; a < 12; ++a) {
        ADDeriv seed = ADDeriv::Zero();
        seed(a) = 1.0;
        u_ad(a) = AD(u_elem(a), seed);
      }

      const ADVec R_loc_ad =
        assemble_element_residual_template<AD>(u_ad, lambda_elem, load);

      for (int a = 0; a < 12; ++a) {
        residual(idx[a]) += R_loc_ad(a).value();

        const ADDeriv& dRa = R_loc_ad(a).derivatives();
        for (int b = 0; b < 12; ++b) {
          const real_t dj = dRa(b);
          if (dj != 0.0) {
            triplets.emplace_back(idx[a], idx[b], dj);
          }
        }
      }
    }

    jacobian.resize(ndof, ndof);
    jacobian.setFromTriplets(triplets.begin(), triplets.end());
    jacobian.makeCompressed();
  }

  void assemble_system(std::vector<std::array<real_t, 3>> load)
  {
    using ADDeriv = Matrix<real_t, 12, 1>;
    using AD = AutoDiffScalar<ADDeriv>;
    using ADVec = Matrix<AD, 12, 1>;
    using Tpl = Triplet<real_t>;

    residual = VectorXd::Zero(ndof);
    std::vector<Tpl> triplets;
    triplets.reserve(elements * 12 * 12);

    for (size_t e = 0; e < elements; ++e) {
      const auto idx = get_element_dof_indices(e);
      const auto lambda_elem = get_element_lambda(e);
      const Matrix<real_t, 12, 1> u_elem = get_element_state(idx);
      const std::array<std::array<real_t, 3>, 2> load_elem = {
        load[e], load[e + 1]
      };

      ADVec u_ad;
      for (int a = 0; a < 12; ++a) {
        ADDeriv seed = ADDeriv::Zero();
        seed(a) = 1.0;
        u_ad(a) = AD(u_elem(a), seed);
      }

      const ADVec R_loc_ad =
        assemble_element_residual_template<AD>(u_ad, lambda_elem, load_elem);

      for (int a = 0; a < 12; ++a) {
        residual(idx[a]) += R_loc_ad(a).value();

        const ADDeriv& dRa = R_loc_ad(a).derivatives();
        for (int b = 0; b < 12; ++b) {
          const real_t dj = dRa(b);
          if (dj != 0.0) {
            triplets.emplace_back(idx[a], idx[b], dj);
          }
        }
      }
    }

    jacobian.resize(ndof, ndof);
    jacobian.setFromTriplets(triplets.begin(), triplets.end());
    jacobian.makeCompressed();
  }

  /**
   * @brief Apply boundary conditions to the residual and jacobian
   *
   * Modifies the residual and jacobian according to boundary conditions:
   * - free_bc: No constraints
   * - simple_bc: Position constraints only
   * - clamped_bc: Position and slope constraints
   */

  void apply_boundary_conditions()
  {
    for (size_t bi = 0; bi < 2; ++bi) {

      EulerBeamBCEnd bcend = boundary_conditions.end[bi];
      size_t ni = 0;
      std::vector<size_t> idx(6);

      switch (bcend) {
        case left:
          ni = 0;
          break;
        case right:
          ni = nodes - 1;
          break;
      }

      EulerBeamBCType bctype = boundary_conditions.type[bi];
      EulerBeamBCVals bcvals = boundary_conditions.vals[bi];
      std::vector<real_t> vals(6);

      switch (bctype) {
        case free_bc:
          idx = {};
          vals = {};
          break;
        case simple_bc:
          idx = { offset_x + 2 * ni + 0,
                  offset_y + 2 * ni + 0,
                  offset_z + 2 * ni + 0 };
          vals = { bcvals.position[0], bcvals.position[1], bcvals.position[2] };
          break;
        case clamped_bc:
          idx = { offset_x + 2 * ni + 0, offset_x + 2 * ni + 1,
                  offset_y + 2 * ni + 0, offset_y + 2 * ni + 1,
                  offset_z + 2 * ni + 0, offset_z + 2 * ni + 1 };
          vals = { bcvals.position[0], bcvals.slope[0],    bcvals.position[1],
                   bcvals.slope[1],    bcvals.position[2], bcvals.slope[2] };
          break;
        case point_force_bc:
          idx = {
            offset_x + 2 * ni + 0,
            offset_y + 2 * ni + 0,
            offset_z + 2 * ni + 0,
          };
          switch (bcend) {
            case left:
              vals = { bcvals.force[0], bcvals.force[1], bcvals.force[2] };
              break;
            case right:
              vals = { bcvals.force[0], bcvals.force[1], bcvals.force[2] };
              break;
          };
          break;
        case point_torque_bc:
          idx = {
            offset_x + 2 * ni + 1,
            offset_y + 2 * ni + 1,
            offset_z + 2 * ni + 1,
          };
          switch (bcend) {
            case left:
              vals = { bcvals.torque[0], bcvals.torque[1], bcvals.torque[2] };
              break;
            case right:
              vals = { bcvals.torque[0], bcvals.torque[1], bcvals.torque[2] };
              break;
          };
          break;
        default:
          break;
      }

      switch (bctype) {
        case point_force_bc:
          for (size_t i = 0; i < vals.size(); i++) {
            // 4) Overwrite the residual to enforce u(i)=vals:
            residual[idx[i]] -= vals[i];
          }
          break;
        case point_torque_bc:
          for (size_t i = 0; i < vals.size(); i++) {
            // 4) Overwrite the residual to enforce u(i)=vals:
            residual[idx[i]] -= vals[i];
          }
          break;
        default:
          for (auto i : idx) {
            // 1) Zero out column i:
            for (SparseMatrix<real_t>::InnerIterator it(jacobian, i); it;
                 ++it) {
              it.valueRef() = 0.0;
            }

            // 2) Zero out row i: because Eigen is column‐major, we loop over
            // each column
            for (int col = 0; col < jacobian.outerSize(); ++col) {
              for (SparseMatrix<real_t>::InnerIterator it(jacobian, col); it;
                   ++it) {
                if (it.row() == i) {
                  it.valueRef() = 0.0;
                }
              }
            }

            // 3) Reinstate the diagonal entry A(i,i) = 1:
            jacobian.coeffRef(i, i) = 1.0;
          }

          for (size_t i = 0; i < vals.size(); i++) {
            // 4) Overwrite the residual to enforce u(i)=vals:
            residual[idx[i]] = u[idx[i]] - vals[i];
          }
          break;
      }
    }
  }

  /**
   * @brief Assemble the residual vector for the nonlinear system
   *
   * @tparam T Scalar type (real_t or autodiff)
   * @param u Current solution vector
   * @param load Current load on the beam
   * @return Assembled residual vector containing:
   *         - Bending energy terms
   *         - External load contributions
   *         - Inextensibility constraints
   *
   * Uses Gauss quadrature with 3 points for numerical integration.
   */
  template<typename T>
  Matrix<T, Dynamic, 1> assemble_residual_template(
    const Matrix<T, Dynamic, 1>& u,
    const std::array<real_t, 3> load) const
  {
    real_t xi_q[] = { 0.1127016654, 0.5, 0.8872983346 };
    real_t w_q[] = { 0.2777777778, 0.4444444444, 0.2777777778 };

    Matrix<T, Dynamic, 1> residual = Matrix<T, Dynamic, 1>::Zero(ndof);

    for (size_t e = 0; e < elements; ++e) {
      std::vector<size_t> elem_nodes = { e, e + 1 };
      std::vector<size_t> idx_x = { offset_x + 2 * elem_nodes[0],
                                    offset_x + 2 * elem_nodes[0] + 1,
                                    offset_x + 2 * elem_nodes[1],
                                    offset_x + 2 * elem_nodes[1] + 1 };
      std::vector<size_t> idx_y = { offset_y + 2 * elem_nodes[0],
                                    offset_y + 2 * elem_nodes[0] + 1,
                                    offset_y + 2 * elem_nodes[1],
                                    offset_y + 2 * elem_nodes[1] + 1 };
      std::vector<size_t> idx_z = { offset_z + 2 * elem_nodes[0],
                                    offset_z + 2 * elem_nodes[0] + 1,
                                    offset_z + 2 * elem_nodes[1],
                                    offset_z + 2 * elem_nodes[1] + 1 };
      std::vector<size_t> idx_l = { elem_nodes[0], elem_nodes[1] };

      std::array<T, 4> ux = {
        u[idx_x[0]], u[idx_x[1]], u[idx_x[2]], u[idx_x[3]]
      };
      std::array<T, 4> uy = {
        u[idx_y[0]], u[idx_y[1]], u[idx_y[2]], u[idx_y[3]]
      };
      std::array<T, 4> uz = {
        u[idx_z[0]], u[idx_z[1]], u[idx_z[2]], u[idx_z[3]]
      };
      std::array<T, 2> ul = { lambda[idx_l[0]], lambda[idx_l[1]] };

      std::vector<T> R_loc_x(4, 0);
      std::vector<T> R_loc_y(4, 0);
      std::vector<T> R_loc_z(4, 0);
      // std::vector<T> R_loc_l(2, 0);

      for (size_t qi = 0; qi < 3; ++qi) {
        real_t xi = xi_q[qi];
        real_t w = w_q[qi];

        auto H = CubicHermite<real_t>::values(xi, ds);
        auto dH = CubicHermite<real_t>::derivs(xi, ds);
        auto ddH = CubicHermite<real_t>::second_derivs(xi, ds);
        auto M = LinearShape<real_t>::values(xi);

        T x = 0, xp = 0, xpp = 0;
        T y = 0, yp = 0, ypp = 0;
        T z = 0, zp = 0, zpp = 0;
        for (size_t i = 0; i < 4; i++) {
          x += H[i] * ux[i];
          xp += dH[i] * ux[i];
          xpp += ddH[i] * ux[i];
          y += H[i] * uy[i];
          yp += dH[i] * uy[i];
          ypp += ddH[i] * uy[i];
          z += H[i] * uz[i];
          zp += dH[i] * uz[i];
          zpp += ddH[i] * uz[i];
        }

        T l = 0;
        for (size_t i = 0; i < 2; i++) {
          l += M[i] * ul[i];
        }

        T S = xp * xp + yp * yp + zp * zp - 1.0;

        for (size_t a = 0; a < 4; ++a) {
          // Bending Energy Contributions
          R_loc_x[a] += EI * xpp * ddH[a] * w * ds;
          R_loc_y[a] += EI * ypp * ddH[a] * w * ds;
          R_loc_z[a] += EI * zpp * ddH[a] * w * ds;
          // External load contribution in y
          R_loc_x[a] -= load[0] * H[a] * w * ds;
          R_loc_y[a] -= load[1] * H[a] * w * ds;
          R_loc_z[a] -= load[2] * H[a] * w * ds;
          // Constraint contributions
          R_loc_x[a] += 2 * (l + r_penalty * S) * xp * dH[a] * w * ds;
          R_loc_y[a] += 2 * (l + r_penalty * S) * yp * dH[a] * w * ds;
          R_loc_z[a] += 2 * (l + r_penalty * S) * zp * dH[a] * w * ds;
        }

        // for (size_t a = 0; a < 2; ++a) {
        //   // Variation w.r.t. lambda
        //   R_loc_l[a] += S * M[a] * w * ds;
        // }
      }

      // Scatter local residual contribution to global residual
      for (size_t i = 0; i < 4; ++i) {
        residual[idx_x[i]] += R_loc_x[i];
        residual[idx_y[i]] += R_loc_y[i];
        residual[idx_z[i]] += R_loc_z[i];
      }

      // for (size_t i = 0; i < 2; ++i) {
      //   residual[idx_l[i]] += R_loc_l[i];
      // }
    }
    return residual;
  }

  template<typename T>
  Matrix<T, Dynamic, 1> assemble_residual_template(
    const Matrix<T, Dynamic, 1>& u,
    const std::vector<std::array<real_t, 3>> load) const
  {
    ELFF_ASSERT(load.size() == nodes,
                "Nodes does not match load vector size.\n");

    real_t xi_q[] = { 0.1127016654, 0.5, 0.8872983346 };
    real_t w_q[] = { 0.2777777778, 0.4444444444, 0.2777777778 };

    Matrix<T, Dynamic, 1> residual = Matrix<T, Dynamic, 1>::Zero(ndof);

    for (size_t e = 0; e < elements; ++e) {
      std::vector<size_t> elem_nodes = { e, e + 1 };
      std::vector<size_t> idx_x = { offset_x + 2 * elem_nodes[0],
                                    offset_x + 2 * elem_nodes[0] + 1,
                                    offset_x + 2 * elem_nodes[1],
                                    offset_x + 2 * elem_nodes[1] + 1 };
      std::vector<size_t> idx_y = { offset_y + 2 * elem_nodes[0],
                                    offset_y + 2 * elem_nodes[0] + 1,
                                    offset_y + 2 * elem_nodes[1],
                                    offset_y + 2 * elem_nodes[1] + 1 };
      std::vector<size_t> idx_z = { offset_z + 2 * elem_nodes[0],
                                    offset_z + 2 * elem_nodes[0] + 1,
                                    offset_z + 2 * elem_nodes[1],
                                    offset_z + 2 * elem_nodes[1] + 1 };
      std::vector<size_t> idx_load = { elem_nodes[0], elem_nodes[1] };
      std::vector<size_t> idx_l = { elem_nodes[0], elem_nodes[1] };

      std::array<T, 4> ux = {
        u[idx_x[0]], u[idx_x[1]], u[idx_x[2]], u[idx_x[3]]
      };
      std::array<T, 4> uy = {
        u[idx_y[0]], u[idx_y[1]], u[idx_y[2]], u[idx_y[3]]
      };
      std::array<T, 4> uz = {
        u[idx_z[0]], u[idx_z[1]], u[idx_z[2]], u[idx_z[3]]
      };
      std::array<real_t, 2> fx = { load[idx_load[0]][0], load[idx_load[1]][0] };
      std::array<real_t, 2> fy = { load[idx_load[0]][1], load[idx_load[1]][1] };
      std::array<real_t, 2> fz = { load[idx_load[0]][2], load[idx_load[1]][2] };
      std::array<T, 2> ul = { lambda[idx_l[0]], lambda[idx_l[1]] };

      std::vector<T> R_loc_x(4, 0);
      std::vector<T> R_loc_y(4, 0);
      std::vector<T> R_loc_z(4, 0);

      for (size_t qi = 0; qi < 3; ++qi) {
        real_t xi = xi_q[qi];
        real_t w = w_q[qi];

        auto H = CubicHermite<real_t>::values(xi, ds);
        auto dH = CubicHermite<real_t>::derivs(xi, ds);
        auto ddH = CubicHermite<real_t>::second_derivs(xi, ds);
        auto M = LinearShape<real_t>::values(xi);

        T x = 0, xp = 0, xpp = 0;
        T y = 0, yp = 0, ypp = 0;
        T z = 0, zp = 0, zpp = 0;
        for (size_t i = 0; i < 4; i++) {
          x += H[i] * ux[i];
          xp += dH[i] * ux[i];
          xpp += ddH[i] * ux[i];
          y += H[i] * uy[i];
          yp += dH[i] * uy[i];
          ypp += ddH[i] * uy[i];
          z += H[i] * uz[i];
          zp += dH[i] * uz[i];
          zpp += ddH[i] * uz[i];
        }

        T l = 0;
        real_t fxp = 0;
        real_t fyp = 0;
        real_t fzp = 0;
        for (size_t i = 0; i < 2; i++) {
          l += M[i] * ul[i];
          fxp += M[i] * fx[i];
          fyp += M[i] * fy[i];
          fzp += M[i] * fz[i];
        }

        T S = xp * xp + yp * yp + zp * zp - 1.0;

        for (size_t a = 0; a < 4; ++a) {
          R_loc_x[a] += EI * xpp * ddH[a] * w * ds;
          R_loc_y[a] += EI * ypp * ddH[a] * w * ds;
          R_loc_z[a] += EI * zpp * ddH[a] * w * ds;
          R_loc_x[a] -= fxp * H[a] * w * ds;
          R_loc_y[a] -= fyp * H[a] * w * ds;
          R_loc_z[a] -= fzp * H[a] * w * ds;
          R_loc_x[a] += 2 * (l + r_penalty * S) * xp * dH[a] * w * ds;
          R_loc_y[a] += 2 * (l + r_penalty * S) * yp * dH[a] * w * ds;
          R_loc_z[a] += 2 * (l + r_penalty * S) * zp * dH[a] * w * ds;
        }
      }

      for (size_t i = 0; i < 4; ++i) {
        residual[idx_x[i]] += R_loc_x[i];
        residual[idx_y[i]] += R_loc_y[i];
        residual[idx_z[i]] += R_loc_z[i];
      }
    }
    return residual;
  }
};

} // namespace ELFF
} // namespace Models
