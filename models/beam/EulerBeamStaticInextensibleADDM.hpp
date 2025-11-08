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
#include <unsupported/Eigen/AutoDiff>

namespace beam {

/**
 * @brief A class to solve the static inextensible Euler–Bernoulli beam
 * equation; this model is valid for large deflections. The strong form of the
 * system that we seek to solve is
 * \f[
 *    \frac{\partial^2}{\partial s^2} \left(EI \frac{\partial^2
 * \mathbf{r}}{\partial s^2}\right) = \mathbf{q}(s)
 * \f]
 * where \f(EI\f) is a bending stiffness, \f(\mathbf{r}(s) = (x(s),y(s))\f) is
 * the deflection of our beam. Unlike the classic Euler-Bernouli beam equation,
 * where inextensibility is implicitly enforced, in the \f(n > 1\f) dimensional
 * version, we consider the deflection of the beam in each dimension along a
 * curvilinear coordinate system, and our inextensbility must be enforced
 * explicitly. Then, in addition, we enforce on the solution a pointwise
 * constraint
 * \f[
 *    ||\mathbf{r}'(s)||^2 = 1.
 * \f]
 * To derive a weak form, we introduce a smooth test function
 * ...
 * \f[
 *    \mathcal{L}_R(x,y,p,q,\lambda_x,\mu) = J(x,y) + \int_0^L
 * \left[\lambda_x(p-x')
 * + \mu(q-y')\right] ds + \frac{r}{2} \int_0^L \left[(p-x')^2 + (q-y')^2\right]
 * ds.
 * \f]
 */

class EulerBeamStaticInextensibleADDM : public EulerBeam
{

public:
  EulerBeamStaticInextensibleADDM(real_t length,
                                  real_t EI,
                                  size_t nodes,
                                  EulerBeamBCs bcs,
                                  real_t r_penalty)
    : EulerBeam(length, EI, nodes, bcs)

    , elements(nodes - 1)
    , dimension(3)
    , dof((2 * nodes))

    , A(Eigen::MatrixXd::Zero(dof, dof))
    , f_x(Eigen::VectorXd::Zero(dof))
    , f_y(Eigen::VectorXd::Zero(dof))
    , f_z(Eigen::VectorXd::Zero(dof))
    , x(Eigen::VectorXd::Zero(dof))
    , y(Eigen::VectorXd::Zero(dof))
    , z(Eigen::VectorXd::Zero(dof))
    , llt()

    , xp(Eigen::VectorXd::Zero(nodes + elements))
    , yp(Eigen::VectorXd::Zero(nodes + elements))
    , zp(Eigen::VectorXd::Zero(nodes + elements))
    , lambda_x(Eigen::VectorXd::Zero(nodes + elements))
    , lambda_y(Eigen::VectorXd::Zero(nodes + elements))
    , lambda_z(Eigen::VectorXd::Zero(nodes + elements))
    , p(Eigen::VectorXd::Ones(nodes + elements))
    , q(Eigen::VectorXd::Zero(nodes + elements))
    , r(Eigen::VectorXd::Zero(nodes + elements))

    , r_penalty(r_penalty)
    , tol_outer(1e-7)
    , max_outer(100000)
  {
    apply_initial_condition();
    assemble_A();
    apply_boundary_condition_A();
    decompose_A();
  };

  ~EulerBeamStaticInextensibleADDM() {};

  /**
   * 
   */
  virtual void solve() {
    solve({0.,0.,0.});
  }

  /**
   * 
   */
  virtual void solve(std::array<real_t,3> load)
  {
    for (size_t iter = 0; iter < max_outer; ++iter) {
      std::cout << "Iter " << iter << std::endl;
      update_pq();
      update_xy(load);
      update_multipliers();
      if (is_converged()) {
        break;
      }
    }

    update_mesh();
  };

  virtual void apply_initial_condition()
  {
    apply_initial_condition_xy();
    apply_initial_condition_pq();
  }

  virtual void apply_initial_condition(EulerBeamMesh& bmesh) override
  {
    apply_initial_condition_xy(bmesh);
    apply_initial_condition_pq();
  }

  const Eigen::VectorXd& get_lambda_x() const { return lambda_x; }
  const Eigen::VectorXd& get_lambda_y() const { return lambda_y; }
  const Eigen::VectorXd& get_lambda_z() const { return lambda_z; }
  const Eigen::MatrixXd& get_A() const { return A; }

protected:
  size_t dimension;
  size_t elements;
  size_t dof;
  real_t r_penalty, alpha;

  Eigen::MatrixXd A;
  Eigen::VectorXd x, y, z;
  Eigen::VectorXd f_x, f_y, f_z;

  Eigen::LLT<Eigen::MatrixXd> llt;

  Eigen::VectorXd lambda_x, lambda_y, lambda_z;
  Eigen::VectorXd p, q, r;
  Eigen::VectorXd xp, yp, zp;

  size_t max_outer;
  real_t tol_outer;

  EulerBeamStaticInextensibleADDM(real_t length,
                                  real_t EI,
                                  real_t mu,
                                  size_t nodes,
                                  EulerBeamBCs bcs,
                                  real_t r_penalty)
    : EulerBeam(length, EI, mu, nodes, bcs)
    , elements(nodes - 1)
    , dimension(3)
    , dof((2 * nodes))

    , A(Eigen::MatrixXd::Zero(dof, dof))
    , f_x(Eigen::VectorXd::Zero(dof))
    , f_y(Eigen::VectorXd::Zero(dof))
    , f_z(Eigen::VectorXd::Zero(dof))
    , x(Eigen::VectorXd::Zero(dof))
    , y(Eigen::VectorXd::Zero(dof))
    , z(Eigen::VectorXd::Zero(dof))
    , llt()

    , xp(Eigen::VectorXd::Zero(nodes + elements))
    , yp(Eigen::VectorXd::Zero(nodes + elements))
    , zp(Eigen::VectorXd::Zero(nodes + elements))
    , lambda_x(Eigen::VectorXd::Zero(nodes + elements))
    , lambda_y(Eigen::VectorXd::Zero(nodes + elements))
    , lambda_z(Eigen::VectorXd::Zero(nodes + elements))
    , p(Eigen::VectorXd::Ones(nodes + elements))
    , q(Eigen::VectorXd::Zero(nodes + elements))
    , r(Eigen::VectorXd::Zero(nodes + elements))

    , r_penalty(r_penalty)
    , tol_outer(1e-7)
    , max_outer(100000)
  {
    apply_initial_condition(mesh);
    assemble_A();
    apply_boundary_condition_A();
    decompose_A();
  };

  // Update this->mesh object based on the solution
  void update_mesh()
  {
    size_t nodes = this->mesh.get_nodes();
    std::vector<std::array<real_t, 3>>& centerline = mesh.get_centerline();
    std::vector<real_t>& s = mesh.get_curvilinear_axis();

    for (size_t i = 0; i < nodes; ++i) {
      centerline[i][0] = x(2 * i);
      centerline[i][1] = y(2 * i);
      centerline[i][2] = z(2 * i);
    }
  }


  // Set initial condition
  void apply_initial_condition_xy()
  {
    real_t h = mesh.get_ds();
    size_t nodes = mesh.get_nodes();

    for (size_t i = 0; i < nodes; i++) {
      x(2 * i + 0) = h * i;
      x(2 * i + 1) = 1.;
      y(2 * i + 0) = 0.;
      y(2 * i + 1) = 0.;
      z(2 * i + 0) = 0.;
      z(2 * i + 1) = 0.;
    }

    update_mesh();
  }

  // Set initial condition
  void apply_initial_condition_xy(EulerBeamMesh& bmesh)
  {
    size_t nodes = bmesh.get_nodes();
    BEAM_ASSERT(nodes == mesh.get_nodes(),
                "Node count of the mesh must match current mesh.\n");
    auto centerline = bmesh.get_centerline();
    auto slopes = bmesh.get_slope();

    for (size_t i = 0; i < nodes; i++) {
      x(2 * i + 0) = centerline[i][0];
      x(2 * i + 1) = slopes[i][0];
      y(2 * i + 0) = centerline[i][1];
      y(2 * i + 1) = slopes[i][1];
      z(2 * i + 0) = centerline[i][2];
      z(2 * i + 1) = slopes[i][2];
    }

    update_mesh();
  }

  //
  void compute_slopes_collocation()
  {
    xp.setZero();
    yp.setZero();
    zp.setZero();
    const real_t h = this->mesh.get_ds();
    const size_t nodes = this->mesh.get_nodes();

    // Node collocation: Take the slope directly from the xy-solution
    for (size_t ni = 0; ni < nodes; ni++) {
      xp[ni] = x[2 * ni + 1];
      yp[ni] = y[2 * ni + 1];
      zp[ni] = z[2 * ni + 1];
    }

    // Midpoint colloation: Evaulate slope via Hermite derivatives at xi = 0.5
    for (size_t ei = 0; ei < elements; ei++) {
      std::array<real_t, 4> dH = CubicHermite<real_t>::derivs(0.5, h);

      size_t edofs[4] = {
        2 * (ei + 0) + 0, 2 * (ei + 0) + 1, 2 * (ei + 1) + 0, 2 * (ei + 1) + 1
      };

      xp[nodes + ei] = dH[0] * x[edofs[0]] + dH[1] * x[edofs[1]] +
                       dH[2] * x[edofs[2]] + dH[3] * x[edofs[3]];
      yp[nodes + ei] = dH[0] * y[edofs[0]] + dH[1] * y[edofs[1]] +
                       dH[2] * y[edofs[2]] + dH[3] * y[edofs[3]];
      zp[nodes + ei] = dH[0] * z[edofs[0]] + dH[1] * z[edofs[1]] +
                       dH[2] * z[edofs[2]] + dH[3] * z[edofs[3]];
    }
  }

  void apply_initial_condition_pq()
  {
    const size_t nodes = this->mesh.get_nodes();
    compute_slopes_collocation();
    for (size_t ci = 0; ci < nodes + elements; ci++) {

      switch (dimension) {
        case 2: {
          real_t xprime = xp[ci];
          real_t yprime = yp[ci];
          real_t den = std::max(1e-14, sqrt(xprime * xprime + yprime * yprime));
          p[ci] = xprime / den;
          q[ci] = yprime / den;
        } break;
        case 3: {
          real_t xprime = xp[ci];
          real_t yprime = yp[ci];
          real_t zprime = zp[ci];
          real_t den = std::max(
            1e-14, sqrt(xprime * xprime + yprime * yprime + zprime * zprime));
          p[ci] = xprime / den;
          q[ci] = yprime / den;
          r[ci] = zprime / den;
        } break;
      }
    }
    apply_boundary_condition_pq();
  }

  // Apply the boundary condition to the p,q variables
  void apply_boundary_condition_pq()
  {
    const size_t nodes = mesh.get_nodes();
    for (size_t bi = 0; bi < 2; bi++) {
      if (boundary_conditions.type[bi] == clamped_bc) {
        switch (boundary_conditions.end[bi]) {
          case left:
            p(0) = boundary_conditions.vals[bi].slope[0];
            q(0) = boundary_conditions.vals[bi].slope[1];
            r(0) = boundary_conditions.vals[bi].slope[2];
            break;
          case right:
            p(nodes - 1) = boundary_conditions.vals[bi].slope[0];
            q(nodes - 1) = boundary_conditions.vals[bi].slope[1];
            r(nodes - 1) = boundary_conditions.vals[bi].slope[2];
            break;
        }
      }
    }
  }

  // Perform a pointwise projection of (p^n, q^n) onto the unit circle to obtain
  // (p^{n+1}, q^{n+1})
  virtual void update_pq()
  {
    const real_t h = mesh.get_ds();
    const size_t nodes = mesh.get_nodes();
    compute_slopes_collocation();

    // #pragma omp parallel for schedule(static)
    for (size_t ci = 0; ci < nodes + elements; ci++) {
      switch (dimension) {
        case 2: {
          const real_t xprime = xp[ci];
          const real_t yprime = yp[ci];
          real_t X = xprime - lambda_x[ci] / r_penalty;
          real_t Y = yprime - lambda_y[ci] / r_penalty;
          real_t norm = std::max(1e-6, sqrt(X * X + Y * Y));
          p[ci] = X / norm;
          q[ci] = Y / norm;
        } break;
        case 3: {
          const real_t xprime = xp[ci];
          const real_t yprime = yp[ci];
          const real_t zprime = zp[ci];
          real_t X = xprime - lambda_x[ci] / r_penalty;
          real_t Y = yprime - lambda_y[ci] / r_penalty;
          real_t Z = zprime - lambda_z[ci] / r_penalty;
          real_t norm = std::max(1e-6, sqrt(X * X + Y * Y + Z * Z));
          p[ci] = X / norm;
          q[ci] = Y / norm;
          r[ci] = Z / norm;
        } break;
      }
    }

    apply_boundary_condition_pq();
  }

  // Apply boundary conditions to the matrix A
  void apply_boundary_condition_A()
  {
    size_t nodes = mesh.get_nodes();
    real_t h = mesh.get_ds();

    for (size_t bi = 0; bi < 2; ++bi) {
      EulerBeamBCEnd bcend = boundary_conditions.end[bi];
      size_t ni = 0;

      std::vector<size_t> idx;   // idx.reserve(2);   // (<=2);
      std::vector<real_t> xvals; // xvals.reserve(2); // (<=2);
      std::vector<real_t> yvals; // yvals.reserve(2); // (<=2);
      std::vector<real_t> zvals; // yvals.reserve(2); // (<=2);

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

      switch (bctype) {
        case free_bc:
          idx = {};
          xvals = {};
          yvals = {};
          zvals = {};
          break;
        case simple_bc:
          idx = { 2 * ni + 0 };
          xvals = { bcvals.position[0] };
          yvals = { bcvals.position[1] };
          zvals = { bcvals.position[2] };
          break;
        case clamped_bc:
          idx = { 2 * ni + 0 };
          xvals = { bcvals.position[0] };
          yvals = { bcvals.position[1] };
          zvals = { bcvals.position[2] };
          break;
        case point_force_bc:
          // idx = {2 * ni + 0};
          // xvals = {bcvals.force[0]};
          // yvals = {bcvals.force[1]};
          // zvals = {bcvals.force[2]};
          idx = {};
          xvals = {};
          yvals = {};
          zvals = {};
          break;
        case point_torque_bc:
          // idx = {2 * ni + 0};
          // xvals = {bcvals.force[0]};
          // yvals = {bcvals.force[1]};
          // zvals = {bcvals.force[2]};
          idx = {};
          xvals = {};
          yvals = {};
          zvals = {};
          break;
      }

      // Set boundary conditions in the matrix A
      // switch(bctype) {
      //   default:
      for (size_t i = 0; i < xvals.size(); i++) {
        A.row(idx[i]).setZero();
        A.col(idx[i]).setZero();
        A(idx[i], idx[i]) = 1.;
      }
      // break;
      // }
    }
  }

  // Assemble the matrix A
  void assemble_A()
  {
    const real_t h = this->mesh.get_ds();
    const real_t xi_q[3] = { 0.1127016654, 0.5, 0.8872983346 };
    const real_t w_q[3] = { 0.2777777778, 0.4444444444, 0.2777777778 };

    // Local stiffness matrices
    Eigen::MatrixXd K4 = Eigen::MatrixXd::Zero(this->dof, this->dof);
    Eigen::MatrixXd K2 = Eigen::MatrixXd::Zero(this->dof, this->dof);

    for (size_t e = 0; e < this->elements; ++e) {
      const size_t edofs[4] = {
        2 * (e + 0) + 0, 2 * (e + 0) + 1, 2 * (e + 1) + 0, 2 * (e + 1) + 1
      };

      real_t K4e[4][4] = { { 0 } };
      real_t K2e[4][4] = { { 0 } };

      // Quadtrature
      for (size_t q = 0; q < 3; ++q) {
        const real_t xi = xi_q[q];
        const real_t w = w_q[q];

        auto H = CubicHermite<real_t>::values(xi, h);
        auto dH = CubicHermite<real_t>::derivs(xi, h);
        auto ddH = CubicHermite<real_t>::second_derivs(xi, h);

        // Accumulate local K4e and K2e
        for (size_t a = 0; a < 4; ++a) {
          for (size_t b = 0; b < 4; ++b) {
            K4e[a][b] += ddH[a] * ddH[b] * w * h;
            K2e[a][b] += dH[a] * dH[b] * w * h;
          }
        }
      }

      // Scatter to global triplets
      for (size_t a = 0; a < 4; ++a) {
        for (size_t b = 0; b < 4; ++b) {
          K4(edofs[a], edofs[b]) += K4e[a][b];
          K2(edofs[a], edofs[b]) += K2e[a][b];
        }
      }
    }

    this->A = this->EI * K4 + this->r_penalty * K2;
  }

  // Cholesky decompose the matrix A
  void decompose_A() { llt.compute(A); }

  // Assemble the rhs vector for x-axis
  void assemble_f(std::array<real_t, 3> load)
  {
    f_x.setZero();
    f_y.setZero();
    f_z.setZero();

    const real_t h = this->mesh.get_ds();
    const size_t nodes = this->mesh.get_nodes();

    const real_t xi_q[3] = { 0.1127016654, 0.5, 0.8872983346 };
    const real_t w_q[3] = { 0.2777777778, 0.4444444444, 0.2777777778 };

    for (size_t e = 0; e < elements; ++e) {
      size_t edofs[4] = {
        2 * (e + 0) + 0, 2 * (e + 0) + 1, 2 * (e + 1) + 0, 2 * (e + 1) + 1
      };
      real_t fxe[4] = { 0, 0, 0, 0 };
      real_t fye[4] = { 0, 0, 0, 0 };
      real_t fze[4] = { 0, 0, 0, 0 };

      for (size_t qi = 0; qi < 3; ++qi) {
        real_t xi = xi_q[qi];
        real_t w = w_q[qi];

        auto L = QuadraticLagrange<real_t>::values(xi);
        auto H = CubicHermite<real_t>::values(xi, h);
        auto dH = CubicHermite<real_t>::derivs(xi, h);

        size_t li = e;
        size_t mi = nodes + e;
        size_t ri = e + 1;

        real_t p_val = L[0] * p[li] + L[1] * p[mi] + L[2] * p[ri];
        real_t lambda_x_val =
          L[0] * lambda_x[li] + L[1] * lambda_x[mi] + L[2] * lambda_x[ri];
        real_t q_val = L[0] * q[li] + L[1] * q[mi] + L[2] * q[ri];
        real_t lambda_y_val =
          L[0] * lambda_y[li] + L[1] * lambda_y[mi] + L[2] * lambda_y[ri];
        real_t r_val = L[0] * r[li] + L[1] * r[mi] + L[2] * r[ri];
        real_t lambda_z_val =
          L[0] * lambda_z[li] + L[1] * lambda_z[mi] + L[2] * lambda_z[ri];

        for (size_t a = 0; a < 4; ++a) {
          fxe[a] += (lambda_x_val + r_penalty * p_val) * dH[a] * w * h;
          fye[a] += (lambda_y_val + r_penalty * q_val) * dH[a] * w * h;
          fze[a] += (lambda_z_val + r_penalty * r_val) * dH[a] * w * h;
          fxe[a] += load[0] * H[a] * w * h;
          fye[a] += load[1] * H[a] * w * h;
          fze[a] += load[2] * H[a] * w * h;
        }
      }
      for (size_t a = 0; a < 4; ++a) {
        f_x(edofs[a]) += fxe[a];
        f_y(edofs[a]) += fye[a];
        f_z(edofs[a]) += fze[a];
      }
    }
  }

  // Apply the boundary condition to the f vector(s)
  void apply_boundary_condition_f()
  {
    size_t nodes = mesh.get_nodes();
    real_t h = mesh.get_ds();

    for (size_t bi = 0; bi < 2; ++bi) {
      EulerBeamBCEnd bcend = boundary_conditions.end[bi];
      size_t ni = 0;

      std::vector<size_t> idx;   // idx.reserve(2);   // (<=2);
      std::vector<real_t> xvals; // xvals.reserve(2); // (<=2);
      std::vector<real_t> yvals; // yvals.reserve(2); // (<=2);
      std::vector<real_t> zvals; // yvals.reserve(2); // (<=2);

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

      switch (bctype) {
        case free_bc:
          idx = {};
          xvals = {};
          yvals = {};
          zvals = {};
          break;
        case simple_bc:
          idx = { 2 * ni + 0 };
          xvals = { bcvals.position[0] };
          yvals = { bcvals.position[1] };
          zvals = { bcvals.position[2] };
          break;
        case clamped_bc:
          idx = { 2 * ni + 0 };           // 2 * ni + 1};
          xvals = { bcvals.position[0] }; //, bcvals.slope[0]};
          yvals = { bcvals.position[1] }; //, bcvals.slope[1]};
          zvals = { bcvals.position[2] }; //, bcvals.slope[1]};
          break;
        case point_force_bc:
          idx = { 2 * ni + 0 };
          xvals = { bcvals.force[0] };
          yvals = { bcvals.force[1] };
          zvals = { bcvals.force[2] };
          break;
        case point_torque_bc:
          idx = { 2 * ni + 0 };
          xvals = { bcvals.torque[0] };
          yvals = { bcvals.torque[1] };
          zvals = { bcvals.torque[2] };
          break;
      }

      // Set boundary conditions in the matrix A
      switch (bctype) {
        case point_force_bc:
          for (size_t i = 0; i < xvals.size(); i++) {
            // Set boundary conditions in f_x
            f_x(idx[i]) += xvals[i];
            // Set boundary conditions in f_y
            f_y(idx[i]) += yvals[i];
            // Set boundary conditions in f_y
            f_z(idx[i]) += zvals[i];
          }
          break;
        default:
          for (size_t i = 0; i < xvals.size(); i++) {
            // Set boundary conditions in f_x
            f_x(idx[i]) = xvals[i];
            // Set boundary conditions in f_y
            f_y(idx[i]) = yvals[i];
            // Set boundary conditions in f_y
            f_z(idx[i]) = zvals[i];
          }
          break;
      }
    }
  }

  // Perform one update of the linear system
  virtual void update_xy(std::array<real_t, 3> load)
  {
    assemble_f(load);
    apply_boundary_condition_f();

    x = llt.solve(f_x);
    y = llt.solve(f_y);
    z = llt.solve(f_z);
  }

  // Update the lagrange multipliers \lambda_x and \mu
  virtual void update_multipliers()
  {
    const real_t h = mesh.get_ds();
    const size_t nodes = mesh.get_nodes();
    compute_slopes_collocation();

    // #pragma omp parallel for schedule(static)
    for (size_t ci = 0; ci < nodes + elements; ci++) {
      lambda_x[ci] += r_penalty * (p[ci] - xp[ci]);
      lambda_y[ci] += r_penalty * (q[ci] - yp[ci]);
      lambda_z[ci] += r_penalty * (r[ci] - zp[ci]);
    }
  }


  // compute the inf norm ||\mathbf{p} - \mathbf{x}'||_{\infty} and ||\mathbf{q}
  // - \mathbf{y}'||_{\infty}
  bool is_converged(bool recompute_slopes = true)
  {
    if (recompute_slopes) {
      compute_slopes_collocation(); // fills xp, yp from current x,y
    }
    const real_t res_p = (p - xp).cwiseAbs().maxCoeff();
    const real_t res_q = (q - yp).cwiseAbs().maxCoeff();
    const real_t res_r = (r - zp).cwiseAbs().maxCoeff();
    const real_t res = std::max({ res_p, res_q, res_r });

    std::cout << "max|p-x'|=" << res_p << ", max|q-y'|=" << res_q << '\n';
    return res < tol_outer;
  }
};

} // namespace beam