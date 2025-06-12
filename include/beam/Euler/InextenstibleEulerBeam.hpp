#pragma once

#include <beam/Euler/EulerBeam.hpp>
#include <beam/FiniteElement/Shapes.hpp>
#include <beam/LinAlg/Matrix.hpp>
#include <beam/LinAlg/Vector.hpp>

// #include <cblas.h>
#include <iostream>
#include <stdlib.h>
#include <vector>

#include <Eigen/Dense>
#include <unsupported/Eigen/AutoDiff>

namespace beam {

// Standalone function to assemble the residual vector
template <class T>
inline Eigen::Matrix<T, Eigen::Dynamic, 1>
_InextensibleEulerBeam_assemble_residual(
    double EI, double load, size_t nodes, double r_penalty, double h,
    const Eigen::Matrix<T, Eigen::Dynamic, 1> &u) {
  double xi_q[] = {0.1127016654, 0.5, 0.8872983346};
  double w_q[] = {0.2777777778, 0.4444444444, 0.2777777778};

  size_t ndof_x = 2 * nodes;
  size_t ndof_y = 2 * nodes;
  size_t ndof_l = nodes;

  size_t offset_x = 0;
  size_t offset_y = ndof_x;
  size_t offset_l = ndof_x + ndof_y;

  size_t elements = nodes - 1;
  size_t dof = ndof_x + ndof_y + ndof_l;

  Eigen::Matrix<T, Eigen::Dynamic, 1> residual =
      Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(dof);

  for (size_t e = 0; e < elements; ++e) {
    std::vector<size_t> elem_nodes = {e, e + 1};
    std::vector<size_t> idx_x = {
        offset_x + 2 * elem_nodes[0], offset_x + 2 * elem_nodes[0] + 1,
        offset_x + 2 * elem_nodes[1], offset_x + 2 * elem_nodes[1] + 1};
    std::vector<size_t> idx_y = {
        offset_y + 2 * elem_nodes[0], offset_y + 2 * elem_nodes[0] + 1,
        offset_y + 2 * elem_nodes[1], offset_y + 2 * elem_nodes[1] + 1};
    std::vector<size_t> idx_l = {offset_l + 1 * elem_nodes[0],
                                 offset_l + 1 * elem_nodes[1]};

    std::array<T, 4> ux = {u[idx_x[0]], u[idx_x[1]], u[idx_x[2]], u[idx_x[3]]};
    std::array<T, 4> uy = {u[idx_y[0]], u[idx_y[1]], u[idx_y[2]], u[idx_y[3]]};
    std::array<T, 2> ul = {u[idx_l[0]], u[idx_l[1]]};

    std::vector<T> R_loc_x(4, 0);
    std::vector<T> R_loc_y(4, 0);
    std::vector<T> R_loc_l(2, 0);

    for (size_t qi = 0; qi < 3; ++qi) {
      double xi = xi_q[qi];
      double w = w_q[qi];

      auto H = CubicHermite<double>::values(xi, h);
      auto dH = CubicHermite<double>::derivs(xi, h);
      auto ddH = CubicHermite<double>::second_derivs(xi, h);
      auto M = LinearShape<double>::values(xi);

      // T x = beam::dot<T,4>(H, ux);
      // T xp = beam::dot<T,4>(dH, ux);
      // T xpp = beam::dot<T,4>(ddH, ux);
      // T y = beam::dot<T,4>(H, uy);
      // T yp = beam::dot<T,4>(dH, uy);
      // T ypp = beam::dot<T,4>(ddH, uy);

      T x = 0, xp = 0, xpp = 0;
      T y = 0, yp = 0, ypp = 0;
      for (size_t i = 0; i < 4; i++) {
        x += H[i] * ux[i];
        xp += dH[i] * ux[i];
        xpp += ddH[i] * ux[i];
        y += H[i] * uy[i];
        yp += dH[i] * uy[i];
        ypp += ddH[i] * uy[i];
      }

      // T l = beam::dot<T,2>(M, ul);
      T l = 0;
      for (size_t i = 0; i < 2; i++) {
        l += M[i] * ul[i];
      }

      T S = xp * xp + yp * yp - 1.0;

      for (size_t a = 0; a < 4; ++a) {
        // Bending Energy Contributions
        R_loc_x[a] += EI * xpp * ddH[a] * w * h;
        R_loc_y[a] += EI * ypp * ddH[a] * w * h;
        // External load contribution in y
        R_loc_y[a] -= load * H[a] * w * h;
        // Constraint contributions
        R_loc_x[a] += 2 * (l + r_penalty * S) * xp * dH[a] * w * h;
        R_loc_y[a] += 2 * (l + r_penalty * S) * yp * dH[a] * w * h;
      }

      for (size_t a = 0; a < 2; ++a) {
        // Variation w.r.t. lambda
        R_loc_l[a] += S * M[a] * w * h;
      }
    }

    // Scatter local residual contribution to global residual
    for (size_t i = 0; i < 4; ++i) {
      residual[idx_x[i]] += R_loc_x[i];
      residual[idx_y[i]] += R_loc_y[i];
    }

    for (size_t i = 0; i < 2; ++i) {
      residual[idx_l[i]] += R_loc_l[i];
    }
  }
  return residual;
}

/// Forward‐difference Jacobian of your residual routine.
/// Calls your double‐precision assembler N times to fill each column.
inline Eigen::MatrixXd _InextensibleEulerBeam_assemble_jacobian_fd(
    double EI, double load, size_t nodes, double r_penalty, double h,
    Eigen::VectorXd const &u, double eps = 1e-6) {
  // number of DOFs
  int dof = int(u.size());
  // pre‐allocate Jacobian
  Eigen::MatrixXd J(dof, dof);
  // compute unperturbed residual
  Eigen::VectorXd R0 = _InextensibleEulerBeam_assemble_residual<double>(
      EI, load, nodes, r_penalty, h, u);
  // finite‐difference each column
  for (int i = 0; i < dof; ++i) {
    // copy & perturb one DOF
    Eigen::VectorXd up = u;
    up(i) += eps;
    // recompute residual
    Eigen::VectorXd R1 = _InextensibleEulerBeam_assemble_residual<double>(
        EI, load, nodes, r_penalty, h, up);
    // fill column i
    J.col(i) = (R1 - R0) / eps;
  }
  return J;
}

class InextensibleEulerBeam {
public:
  InextensibleEulerBeam(double length, double EI, double load, double area,
                        size_t nodes, double r_penalty, beam::EulerBeamBCs bcs)
      : length(length), EI(EI), load(load), area(area), nodes(nodes),
        dimension(2), elements(nodes - 1), dof((2 * nodes * dimension) + nodes),
        h(length / elements), r_penalty(r_penalty), boundary_conditions(bcs),
        jacobian(dof, dof), residual(dof), u(dof) {};

  ~InextensibleEulerBeam() {};

  const Eigen::VectorXd &get_solution() const { return u; }
  const Eigen::VectorXd &get_residual() const { return residual; }
  const Eigen::MatrixXd &get_jacobian() const { return jacobian; }

  void set_solution(Eigen::MatrixXd sol) { this->u = sol; }
  void set_residual(Eigen::MatrixXd r) { this->residual = r; }
  void set_jacobian(Eigen::MatrixXd j) { this->jacobian = j; }

  void assemble_residual() {
    residual = beam::_InextensibleEulerBeam_assemble_residual<double>(
        EI, load, nodes, r_penalty, h, u);
  }

  void assemble_jacobian_fd(double eps = 1e-6) {
    int dof = int(u.size());
    Eigen::VectorXd R0 = _InextensibleEulerBeam_assemble_residual<double>(EI, load, nodes, r_penalty, h, u);
    for (int i = 0; i < dof; ++i) {
      Eigen::VectorXd up = u;
      up(i) += eps;
      Eigen::VectorXd R1 = _InextensibleEulerBeam_assemble_residual<double>(
          EI, load, nodes, r_penalty, h, up);
      jacobian.col(i) = (R1 - R0) / eps;
    }
  }

  void assemble_with_jacobian() {
    using namespace Eigen;
    // 2) Define the AutoDiff scalar & vector types
    using Grad = Eigen::Matrix<double, Eigen::Dynamic, 1>;
    using AD = Eigen::AutoDiffScalar<Grad>;
    using ADVec = Eigen::Matrix<AD, Dynamic, 1>;

    // 3) seed the input
    ADVec x_ad(dof);
    for (int i = 0; i < dof; ++i) {
      Grad seed = Grad::Zero(dof);
      seed(i) = 1.0;
      x_ad(i) = AD(u(i), seed);
    }

    // 4) call the templated residual
    ADVec R_ad = beam::_InextensibleEulerBeam_assemble_residual<AD>(
        EI, load, nodes, r_penalty, h, x_ad);

    // 5) unpack values + derivatives
    residual.resize(dof);
    jacobian.resize(dof, dof);
    for (int i = 0; i < dof; ++i) {
      residual(i) = R_ad(i).value();
      jacobian.row(i) = R_ad(i).derivatives().transpose();
    }
  }

private:
  double length;
  double EI;
  double load;
  double area;
  size_t nodes;
  size_t dimension;
  size_t elements;
  size_t dof;
  double h;
  double r_penalty;
  beam::EulerBeamBCs boundary_conditions;

  // Matrix<double> jacobian;
  Eigen::VectorXd residual;
  Eigen::MatrixXd jacobian;
  Eigen::VectorXd u;
};

} // namespace beam