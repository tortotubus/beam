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

namespace beam {

template <class T, std::size_t Size> T dot(std::array<T,Size> left, std::array<T,Size> right) {
  T val = 0;
  for (size_t i = 0; i < left.size(); i++) {
    val += left[i] * right[i];
  }

  return val;
}

// Standalone function to assemble the residual vector
template <class T>
Eigen::VectorXd _InextensibleEulerBeam_assemble_residual(
    T length, T EI, T load, size_t nodes, T r_penalty, T h,
    const Eigen::VectorXd& u) {
  T xi_q[] = {0.1127016654, 0.5, 0.8872983346};
  T w_q[] = {0.2777777778, 0.4444444444, 0.2777777778};

  size_t ndof_x = 2 * nodes;
  size_t ndof_y = 2 * nodes;
  size_t ndof_l = nodes;

  size_t offset_x = 0;
  size_t offset_y = ndof_x;
  size_t offset_l = ndof_x + ndof_y;

  size_t elements = nodes - 1;
  size_t dof = ndof_x + ndof_y + ndof_l;
  Eigen::VectorXd residual = Eigen::VectorXd::Zero(dof);

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
      T xi = xi_q[qi];
      T w = w_q[qi];

      auto H = CubicHermite<T>::values(xi, h);
      auto dH = CubicHermite<T>::derivs(xi, h);
      auto ddH = CubicHermite<T>::second_derivs(xi, h);
      auto M = LinearShape<T>::values(xi);

      T x = beam::dot<T,4>(H, ux);
      T xp = beam::dot<T,4>(dH, ux);
      T xpp = beam::dot<T,4>(ddH, ux);
      T y = beam::dot<T,4>(H, uy);
      T yp = beam::dot<T,4>(dH, uy);
      T ypp = beam::dot<T,4>(ddH, uy);
      T l = beam::dot<T,2>(M, ul);

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



class InextensibleEulerBeam {
public:
  InextensibleEulerBeam(double length, double EI, double load, double area,
                        size_t nodes, double r_penalty, beam::EulerBeamBCs bcs)
      : length(length), EI(EI), load(load), area(area), nodes(nodes),
        dimension(2), elements(nodes - 1), dof((2 * nodes * dimension) + nodes),
        h(length / elements), r_penalty(r_penalty), boundary_conditions(bcs),
        // jacobian(dof, dof, ColMajorOrder),
        jacobian(dof, dof),
        residual(dof), u(dof) {};

  ~InextensibleEulerBeam() {};

  Eigen::VectorXd get_solution() { return u; }
  Eigen::VectorXd get_residual() { return residual; }
  // Matrix<double> get_jacobian() { return jacobian; }

  void set_solution(Eigen::VectorXd sol) { this->u = sol; }
  void set_residual(Eigen::VectorXd r) { this->residual = r; }

  void assemble_residual() {
    residual = beam::_InextensibleEulerBeam_assemble_residual<double>(length, EI, load, nodes, r_penalty, h, u);
  }

  // void assemble_residual() {
  //   double xi_q[] = {0.1127016654, 0.5, 0.8872983346};
  //   double w_q[] = {0.2777777778, 0.4444444444, 0.2777777778};

  //   size_t ndof_x = 2 * nodes;
  //   size_t ndof_y = 2 * nodes;
  //   size_t ndof_l = nodes;

  //   size_t offset_x = 0;
  //   size_t offset_y = ndof_x;
  //   size_t offset_l = ndof_x + ndof_y;

  //   for (size_t e = 0; e < this->elements; ++e) {
  //     std::vector<size_t> nodes = {e, e + 1};
  //     std::vector<size_t> idx_x = {
  //         offset_x + 2 * nodes[0], offset_x + 2 * nodes[0] + 1,
  //         offset_x + 2 * nodes[1], offset_x + 2 * nodes[1] + 1};
  //     std::vector<size_t> idx_y = {
  //         offset_y + 2 * nodes[0], offset_y + 2 * nodes[0] + 1,
  //         offset_y + 2 * nodes[1], offset_y + 2 * nodes[1] + 1};
  //     std::vector<size_t> idx_l = {offset_l + 1 * nodes[0],
  //                                  offset_l + 1 * nodes[1]};

  //     std::array<double, 4> ux = {u[idx_x[0]], u[idx_x[1]], u[idx_x[2]],
  //                                 u[idx_x[3]]};
  //     std::array<double, 4> uy = {u[idx_y[0]], u[idx_y[1]], u[idx_y[2]],
  //                                 u[idx_y[3]]};
  //     std::array<double, 2> ul = {u[idx_l[0]], u[idx_l[1]]};

  //     std::vector<double> R_loc_x(4);
  //     std::vector<double> R_loc_y(4);
  //     std::vector<double> R_loc_l(2);

  //     for (size_t qi = 0; qi < 3; ++qi) {
  //       double xi = xi_q[qi];
  //       double w = w_q[qi];
  //       double h = this->h;

  //       auto H = CubicHermite::values(xi, h);
  //       auto dH = CubicHermite::derivs(xi, h);
  //       auto ddH = CubicHermite::second_derivs(xi, h);
  //       auto M = LinearShape::values(xi);

  //       // double x = cblas_ddot(4, H.data(), 1, ux.data(), 1);
  //       // double xp = cblas_ddot(4, dH.data(), 1, ux.data(), 1);
  //       // double xpp = cblas_ddot(4, ddH.data(), 1, ux.data(), 1);

  //       // double y = cblas_ddot(4, H.data(), 1, uy.data(), 1);
  //       // double yp = cblas_ddot(4, dH.data(), 1, uy.data(), 1);
  //       // double ypp = cblas_ddot(4, ddH.data(), 1, uy.data(), 1);

  //       // double l = cblas_ddot(2, M.data(), 1, ul.data(), 1);

  //       double S = xp * xp + yp * yp - 1.0;

  //       for (size_t a = 0; a < 4; ++a) {
  //         // Bending Energy Contributions
  //         R_loc_x[a] += EI * xpp * ddH[a] * w * h;
  //         R_loc_y[a] += EI * ypp * ddH[a] * w * h;
  //         // External load contribution in y
  //         R_loc_y[a] -= load * H[a] * w * h;
  //         // Constraint contributions
  //         R_loc_x[a] += 2 * (l + r_penalty * S) * xp * dH[a] * w * h;
  //         R_loc_y[a] += 2 * (l + r_penalty * S) * yp * dH[a] * w * h;
  //       }

  //       for (size_t a = 0; a < 2; ++a) {
  //         // Variation w.r.t. lambda
  //         R_loc_l[a] += S * M[a] * w * h;
  //       }
  //     }

  //     // Scatter local residual contribution to global residual
  //     for (size_t i = 0; i < 4; ++i) {
  //       residual[idx_x[i]] += R_loc_x[i];
  //       residual[idx_y[i]] += R_loc_y[i];
  //     }

  //     for (size_t i = 0; i < 2; ++i) {
  //       residual[idx_l[i]] += R_loc_l[i];
  //     }
  //   }
  // }

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