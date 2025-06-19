#pragma once


//#include <beam/LinAlg/Matrix.hpp>
//#include <beam/LinAlg/Vector.hpp>

#include "EulerBeam.hpp"

#include <cstdio> // for popen, pclose, fprintf
#include <iostream>
#include <stdlib.h>
#include <utility>
#include <vector>

#include <Eigen/Dense>
#include <unsupported/Eigen/AutoDiff>

namespace beam {

class EulerBeamStaticInextensible {
public:
  EulerBeamStaticInextensible(double length, double EI, double load, double area,
                        size_t nodes, double r_penalty, beam::EulerBeamBCs bcs)
      : length(length), EI(EI), load(load), area(area), nodes(nodes),
        dimension(2), elements(nodes - 1), dof((2 * nodes * dimension) + nodes),
        h(length / elements), r_penalty(r_penalty), boundary_conditions(bcs),
        jacobian(Eigen::MatrixXd::Zero(dof, dof)),
        residual(Eigen::VectorXd::Zero(dof)), u(Eigen::VectorXd::Zero(dof)) {};

  ~EulerBeamStaticInextensible() {};

  const Eigen::VectorXd &get_solution() const { return u; }
  const Eigen::VectorXd &get_residual() const { return residual; }
  const Eigen::MatrixXd &get_jacobian() const { return jacobian; }

  void set_solution(Eigen::MatrixXd sol) { this->u = sol; }
  void set_residual(Eigen::MatrixXd r) { this->residual = r; }
  void set_jacobian(Eigen::MatrixXd j) { this->jacobian = j; }

  void solve() {
    apply_initial_condition();
    double tol = 1e-6;
    size_t max_iter = 100;

    Eigen::LLT<Eigen::MatrixXd> cholesky;
    Eigen::LDLT<Eigen::MatrixXd> ldlt;

    for (size_t it = 0; it < max_iter; it++) {
      assemble_system();
      apply_boundary_conditions();

      double res_norm = residual.norm();
      std::cout << "iter " << it << "  ‖R‖ = " << res_norm << "\n";
      if (res_norm < tol) {
        std::cout << "Converged in " << it << " iters.\n";
        break;
      }

      Eigen::VectorXd delta_u;
      ldlt.compute(jacobian);
      if (1 == 1) { // ldlt.info() == Eigen::Success) {
        delta_u = jacobian.colPivHouseholderQr().solve(-residual);
      } else {
        delta_u = ldlt.solve(-residual);
      }

      // 3) solve J Δu = −R
      // Option A: PartialPivLU (fast, pivoting)
      // Eigen::VectorXd delta_u = jacobian.partialPivLu().solve(-residual);

      // Option B: ColPivHouseholderQR (more robust for rank-deficient)
      // Eigen::VectorXd delta_u =
      // jacobian.colPivHouseholderQr().solve(-residual);

      // 3) solve J Δu = −R
      // Eigen::ConjugateGradient<Eigen::MatrixXd> cg;
      // cg.setMaxIterations(1000);
      // cg.setTolerance(1e-6);
      // Eigen::VectorXd delta_u = cg.compute(jacobian).solve(-residual);
      //
      // if(cg.info() != Eigen::Success) {
      //    std::cout << "CG solver failed to converge in " << cg.iterations()
      //              << " iterations. (error: " << cg.error() << ")\n";
      //}

      // 4) update
      u += delta_u;
    }
  }

  Eigen::MatrixXd get_centerline() {
    Eigen::MatrixXd coords = Eigen::MatrixXd::Zero(nodes, dimension);

    size_t ndof_x = 2 * nodes;
    size_t ndof_y = 2 * nodes;

    size_t offset_x = 0;
    size_t offset_y = ndof_x;

    for (size_t i = 0; i < nodes; ++i) {
      coords(i, 0) = u(offset_x + 2 * i);
      coords(i, 1) = u(offset_y + 2 * i);
    }

    return coords;
  }

  void get_centerline_gnuplot() {
    FILE *pipe = popen("gnuplot -persist", "w");
    if (!pipe) {
      std::cerr << "Failed to open pipe to gnuplot" << std::endl;
      return;
    }

    // Get the beam centerline coordinates
    Eigen::MatrixXd coords = get_centerline();

    // Configure the plot
    fprintf(pipe, "set title 'Beam Centerline'\n");
    fprintf(pipe, "set xlabel 'x'\n");
    fprintf(pipe, "set ylabel 'y'\n");
    fprintf(pipe, "set grid\n");
    fprintf(pipe, "set size ratio -1\n"); // Equal aspect ratio
    fprintf(pipe, "plot '-' using 1:2 with linespoints title 'beam' pt 7\n");

    // Send the data points
    for (Eigen::Index i = 0; i < coords.rows(); ++i) {
      fprintf(pipe, "%f %f\n", coords(i, 0), coords(i, 1));
    }
    fprintf(pipe, "e\n"); // End of data marker for gnuplot

    // Clean up
    pclose(pipe);
  }

  void get_jacobian_gnuplot() {
    FILE *pipe = popen("gnuplot -persist", "w");
    if (!pipe) {
      std::cerr << "Failed to open pipe to gnuplot" << std::endl;
      return;
    }

    // Configure the plot
    fprintf(pipe, "set title 'Jacobian Matrix Heatmap'\n");
    fprintf(pipe, "set xlabel 'Column'\n");
    fprintf(pipe, "set ylabel 'Row'\n");
    fprintf(pipe, "set view map\n");
    fprintf(pipe, "set size square\n");
    fprintf(pipe, "set palette defined (-1 'blue', 0 'white', 1 'red')\n");
    fprintf(pipe, "set cbrange [-1:1]\n");
    fprintf(pipe, "plot '-' matrix with image notitle\n");

    // Normalize the matrix for better visualization
    double max_abs_val = jacobian.cwiseAbs().maxCoeff();
    Eigen::MatrixXd normalized = jacobian / max_abs_val;

    // Send the matrix data
    for (Eigen::Index i = 0; i < jacobian.rows(); ++i) {
      for (Eigen::Index j = 0; j < jacobian.cols(); ++j) {
        fprintf(pipe, "%f ", normalized(i, j));
      }
      fprintf(pipe, "\n");
    }
    fprintf(pipe, "e\n"); // End of data marker for gnuplot

    // Clean up
    pclose(pipe);
  }

private:
  void assemble_residual() { residual = assemble_residual_impl<double>(u); }

  void assemble_jacobian_fd(double eps = 1e-6) {
    Eigen::VectorXd R0 = assemble_residual_impl<double>(u);
    for (int i = 0; i < int(dof); ++i) {
      Eigen::VectorXd up = u;
      up(i) += eps;
      Eigen::VectorXd R1 = assemble_residual_impl<double>(up);
      jacobian.col(i) = (R1 - R0) / eps;
    }
  }

  void assemble_system() {
    using AD = Eigen::AutoDiffScalar<Eigen::VectorXd>;
    using ADVec = Eigen::Matrix<AD, Eigen::Dynamic, 1>;

    ADVec x_ad(dof);
    for (int i = 0; i < int(dof); ++i) {
      Eigen::VectorXd seed = Eigen::VectorXd::Zero(dof);
      seed(i) = 1.0;
      x_ad(i) = AD(u(i), seed);
    }

    ADVec R_ad = assemble_residual_impl<AD>(x_ad);

    residual.resize(dof);
    jacobian.resize(dof, dof);
    for (int i = 0; i < int(dof); ++i) {
      residual(i) = R_ad(i).value();
      jacobian.row(i) = R_ad(i).derivatives().transpose();
    }
  }

  void apply_initial_condition() {
    size_t ndof_x = 2 * nodes;
    size_t ndof_y = 2 * nodes;
    size_t ndof_l = nodes;
    size_t offset_x = 0;
    size_t offset_y = ndof_x;
    size_t offset_l = ndof_x + ndof_y;
    for (size_t i = 0; i < nodes; i++) {
      u(offset_x + 2 * i + 0) = h * i;
      u(offset_x + 2 * i + 1) = h;
      u(offset_y + 2 * i + 0) = 0.;
      u(offset_y + 2 * i + 1) = 0.;
      u(offset_l + i) = 0.;
    }
  }

  /**
   * @brief Apply boundary conditions to the residual and jacobian
   * 
   * Modifies the residual and jacobian according to boundary conditions:
   * - free_bc: No constraints
   * - simple_bc: Position constraints only
   * - clamped_bc: Position and slope constraints
   */
  void apply_boundary_conditions() {
    size_t ndof_x = 2 * nodes;
    size_t ndof_y = 2 * nodes;
    size_t ndof_l = nodes;
    size_t offset_x = 0;
    size_t offset_y = ndof_x;
    size_t offset_l = ndof_x + ndof_y;

    for (size_t bi = 0; bi < 2; ++bi) {
      EulerBeamBCEnd bcend = boundary_conditions.end[bi];
      size_t ni = 0;
      std::vector<size_t> idx(5);

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
      std::vector<double> vals(5);

      switch (bctype) {
      case free_bc:
        idx = {};
        vals = {};
        break;
      case simple_bc:
        idx = {offset_x + 2 * ni + 0, offset_y + 2 * ni + 0};
        vals = {bcvals.position[0], bcvals.position[1]};
        break;
      case clamped_bc:
        idx = {offset_x + 2 * ni + 0, offset_x + 2 * ni + 1,
               offset_y + 2 * ni + 0, offset_y + 2 * ni + 1};
        vals = {bcvals.position[0], bcvals.slope[0], bcvals.position[1],
                bcvals.slope[1]};
        break;
      }

      for (size_t i = 0; i < vals.size(); ++i) {
        residual(idx[i]) = u(idx[i]) - vals[i];
        for (size_t j = 0; j < dof; ++j) {
          jacobian(idx[i], j) = 0.;
          jacobian(j, idx[i]) = 0.;
        }
      }
    }
  }

  /**
   * @brief Assemble the residual vector for the nonlinear system
   * 
   * @tparam T Scalar type (double or autodiff)
   * @param u Current solution vector
   * @return Assembled residual vector containing:
   *         - Bending energy terms
   *         - External load contributions
   *         - Inextensibility constraints
   * 
   * Uses Gauss quadrature with 3 points for numerical integration.
   */
  template <typename T>
  Eigen::Matrix<T, Eigen::Dynamic, 1>
  assemble_residual_impl(const Eigen::Matrix<T, Eigen::Dynamic, 1> &u) const {
    double xi_q[] = {0.1127016654, 0.5, 0.8872983346};
    double w_q[] = {0.2777777778, 0.4444444444, 0.2777777778};

    size_t ndof_x = 2 * nodes;
    size_t ndof_y = 2 * nodes;
    size_t ndof_l = nodes;

    size_t offset_x = 0;
    size_t offset_y = ndof_x;
    size_t offset_l = ndof_x + ndof_y;

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

      std::array<T, 4> ux = {u[idx_x[0]], u[idx_x[1]], u[idx_x[2]],
                             u[idx_x[3]]};
      std::array<T, 4> uy = {u[idx_y[0]], u[idx_y[1]], u[idx_y[2]],
                             u[idx_y[3]]};
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