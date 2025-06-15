#pragma once

#include <beam/Euler/EulerBeam.hpp>
#include <beam/FiniteElement/Shapes.hpp>
#include <beam/LinAlg/Matrix.hpp>
#include <beam/LinAlg/Vector.hpp>

// #include <cblas.h>
#include <cstdio> // for popen, pclose, fprintf
#include <iostream>
#include <stdlib.h>
#include <utility>
#include <vector>

#include <Eigen/Dense>
#include <unsupported/Eigen/AutoDiff>

namespace beam {

class InextensibleEulerBeamFaster {
public:
  InextensibleEulerBeamFaster(double length, double EI, double load,
                              double area, size_t nodes, double r_penalty,
                              beam::EulerBeamBCs bcs)
      : length(length), EI(EI), load(load), area(area), nodes(nodes),
        dimension(2), elements(nodes - 1), dof((2 * nodes * dimension) + nodes),
        h(length / elements), r_penalty(r_penalty), boundary_conditions(bcs),
        jacobian(Eigen::MatrixXd::Zero(dof, dof)),
        residual(Eigen::VectorXd::Zero(dof)), u(Eigen::VectorXd::Zero(dof)) {};

  ~InextensibleEulerBeamFaster() {};

  const Eigen::VectorXd &get_solution() const { return u; }
  const Eigen::VectorXd &get_residual() const { return residual; }
  const Eigen::MatrixXd &get_jacobian() const { return jacobian; }

  void set_solution(Eigen::MatrixXd sol) { this->u = sol; }
  void set_residual(Eigen::MatrixXd r) { this->residual = r; }
  void set_jacobian(Eigen::MatrixXd j) { this->jacobian = j; }

  void assemble_residual_xy() {
    residual = assemble_residual_impl_space<double>(u);
  }

  void assemble_system_xy() {
    using AD = Eigen::AutoDiffScalar<Eigen::VectorXd>;
    using ADVec = Eigen::Matrix<AD, Eigen::Dynamic, 1>;

    ADVec x_ad(dof);
    for (int i = 0; i < int(dof); ++i) {
      Eigen::VectorXd seed = Eigen::VectorXd::Zero(dof);
      seed(i) = 1.0;
      x_ad(i) = AD(u(i), seed);
    }

    ADVec R_ad = assemble_residual_impl_space<AD>(x_ad);

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
    for (size_t i = 0; i < nodes; i++) {
      u(offset_x + 2 * i + 0) = h * i;
      u(offset_x + 2 * i + 1) = h;
      u(offset_y + 2 * i + 0) = 0.;
      u(offset_y + 2 * i + 1) = 0.;
    }
  }

  void apply_boundary_conditions() {
    size_t ndof_x = 2 * nodes;
    size_t ndof_y = 2 * nodes;
    // size_t ndof_l = nodes;
    size_t offset_x = 0;
    size_t offset_y = ndof_x;
    // size_t offset_l = ndof_x + ndof_y;

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

  void solve() {
    apply_initial_condition();
    double tol = 1e-6;
    size_t max_iter = 15;

    for (size_t it = 0; it < max_iter; it++) {
      assemble_system_xy();
      apply_boundary_conditions();

      double res_norm = residual.norm();
      std::cout << "iter " << it << "  ‖R‖ = " << res_norm << "\n";
      if (res_norm < tol) {
        std::cout << "Converged in " << it << " iters.\n";
        break;
      }

      Eigen::VectorXd delta_u;
      delta_u = jacobian.colPivHouseholderQr().solve(-residual);

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
  template <typename T>
  Eigen::Matrix<T, Eigen::Dynamic, 1> assemble_residual_impl_space(
      const Eigen::Matrix<T, Eigen::Dynamic, 1> &u) const {
    double xi_q[] = {0.1127016654, 0.5, 0.8872983346};
    double w_q[] = {0.2777777778, 0.4444444444, 0.2777777778};

    size_t ndof_x = 2 * nodes;
    size_t offset_x = 0;
    size_t offset_y = ndof_x;

    Eigen::Matrix<T, Eigen::Dynamic, 1> residual =
        Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(dof);

    for (size_t e = 0; e < elements; ++e) {

      std::vector<size_t> elem_nodes = {e, e + 1};
      std::vector<size_t> idx(dimension * 2);
      std::array<T, 4> uX(dimension*2);

      for (size_t d = 0; d < 4*dimension; ++d) {
        size_t offset_X = d*2;
        size_t stride_X = dimension*2;
        std::array<size_t,4> idx_d = {
            offset_X + stride_X * elem_nodes[0],     // {U_d^e)
            offset_X + stride_X * elem_nodes[0] + 1, // {U_d^e}^\prime
            offset_X + stride_X * elem_nodes[1],     // {U_d^{e+1}}
            offset_X + stride_X * elem_nodes[1] + 1  // {U_d^{e+1}}^\prime
        };

        idx[d*4 + 0] = id_d[0];
        idx[d*4 + 1] = id_d[1];
        idx[d*4 + 2] = id_d[2];
        idx[d*4 + 3] = id_d[3];

        uX[idx_d[0]] = u[idx_d[0]];
        uX[idx_d[1]] = u[idx_d[1]];
        uX[idx_d[2]] = u[idx_d[2]];
        uX[idx_d[3]] = u[idx_d[3]];
      }

      std::vector<T> R_loc(4, 0);

      for (size_t qi = 0; qi < 3; ++qi) {
        double xi = xi_q[qi];
        double w = w_q[qi];

        auto H = CubicHermite<double>::values(xi, h);
        auto dH = CubicHermite<double>::derivs(xi, h);
        auto ddH = CubicHermite<double>::second_derivs(xi, h);

        T X = 0, Xp = 0, Xpp = 0;
        for (size_t i = 0; i < 4; i++) {
          X += H[i] * uX[i];
          Xp += dH[i] * uX[i];
          Xpp += ddH[i] * uX[i];
        }

        for (size_t i = 0; i < 4; ++i) {
          // Bending Energy Contributions
          R_loc[i] += EI * Xpp * ddH[i] * w * h;
        }
      }

      // Scatter local residual contribution to global residual
      for (size_t i = 0; i < 4; ++i) {
        residual[idx[i]] += R_loc[i];
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

  Eigen::VectorXd residual;
  Eigen::MatrixXd jacobian;
  Eigen::VectorXd u;
};

} // namespace beam