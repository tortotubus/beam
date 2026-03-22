#include "models/beam/EulerBeamStaticInextensibleADDM.hpp"

namespace ELFF {
namespace Models {

EulerBeamStaticInextensibleADDM::EulerBeamStaticInextensibleADDM(
  real_t length,
  real_t EI,
  size_t nodes,
  EulerBeamBCs bcs,
  real_t r_penalty)
  : EulerBeam(length, EI, nodes, bcs)
  , elements(nodes - 1)
  , dimension(3)
  , dof((2 * nodes))
  , r_penalty(r_penalty)
  , alpha()
  , A(MatrixXd::Zero(dof, dof))
  , A_unconstrained(MatrixXd::Zero(dof, dof))
  , x(VectorXd::Zero(dof))
  , y(VectorXd::Zero(dof))
  , z(VectorXd::Zero(dof))
  , f_x(VectorXd::Zero(dof))
  , f_y(VectorXd::Zero(dof))
  , f_z(VectorXd::Zero(dof))
  , llt()
  , lambda_x(VectorXd::Zero(nodes + elements))
  , lambda_y(VectorXd::Zero(nodes + elements))
  , lambda_z(VectorXd::Zero(nodes + elements))
  , p(VectorXd::Ones(nodes + elements))
  , q(VectorXd::Zero(nodes + elements))
  , r(VectorXd::Zero(nodes + elements))
  , xp(VectorXd::Zero(nodes + elements))
  , yp(VectorXd::Zero(nodes + elements))
  , zp(VectorXd::Zero(nodes + elements))
  , max_outer(100000)
  , tol_outer(1e-7)
{
  apply_initial_condition();
  assemble_A();
  apply_boundary_condition_A();
  decompose_A();
}

EulerBeamStaticInextensibleADDM::EulerBeamStaticInextensibleADDM(
  real_t length,
  real_t EI,
  real_t mu,
  size_t nodes,
  EulerBeamBCs bcs,
  real_t r_penalty)
  : EulerBeam(length, EI, mu, nodes, bcs)
  , elements(nodes - 1)
  , dimension(3)
  , dof((2 * nodes))
  , r_penalty(r_penalty)
  , alpha()
  , A(MatrixXd::Zero(dof, dof))
  , A_unconstrained(MatrixXd::Zero(dof, dof))
  , x(VectorXd::Zero(dof))
  , y(VectorXd::Zero(dof))
  , z(VectorXd::Zero(dof))
  , f_x(VectorXd::Zero(dof))
  , f_y(VectorXd::Zero(dof))
  , f_z(VectorXd::Zero(dof))
  , llt()
  , lambda_x(VectorXd::Zero(nodes + elements))
  , lambda_y(VectorXd::Zero(nodes + elements))
  , lambda_z(VectorXd::Zero(nodes + elements))
  , p(VectorXd::Ones(nodes + elements))
  , q(VectorXd::Zero(nodes + elements))
  , r(VectorXd::Zero(nodes + elements))
  , xp(VectorXd::Zero(nodes + elements))
  , yp(VectorXd::Zero(nodes + elements))
  , zp(VectorXd::Zero(nodes + elements))
  , max_outer(100000)
  , tol_outer(1e-7)
{
  apply_initial_condition(mesh);
  assemble_A();
  apply_boundary_condition_A();
  decompose_A();
}

EulerBeamStaticInextensibleADDM::~EulerBeamStaticInextensibleADDM() = default;

void
EulerBeamStaticInextensibleADDM::solve()
{
  solve({ 0., 0., 0. });
}

void
EulerBeamStaticInextensibleADDM::solve(std::array<real_t, 3> load)
{
  for (size_t iter = 0; iter < max_outer; ++iter) {
    update_pq();
    update_xy(load);
    update_multipliers();
    if (is_converged()) {
      break;
    }
  }

  update_mesh();
}

void
EulerBeamStaticInextensibleADDM::apply_initial_condition()
{
  apply_initial_condition_xy();
  apply_initial_condition_pq();
}

void
EulerBeamStaticInextensibleADDM::apply_initial_condition(EulerBeamMesh& bmesh)
{
  apply_initial_condition_xy(bmesh);
  apply_initial_condition_pq();
}

const VectorXd&
EulerBeamStaticInextensibleADDM::get_lambda_x() const
{
  return lambda_x;
}

const VectorXd&
EulerBeamStaticInextensibleADDM::get_lambda_y() const
{
  return lambda_y;
}

const VectorXd&
EulerBeamStaticInextensibleADDM::get_lambda_z() const
{
  return lambda_z;
}

const MatrixXd&
EulerBeamStaticInextensibleADDM::get_A() const
{
  return A;
}

void
EulerBeamStaticInextensibleADDM::update_mesh()
{
  const size_t nodes = this->mesh.get_nodes();
  std::vector<std::array<real_t, 3>>& centerline = mesh.get_centerline();
  std::vector<std::array<real_t, 3>>& slope = mesh.get_slope();
  std::vector<real_t>& s = mesh.get_curvilinear_axis();
  (void)s;

  for (size_t i = 0; i < nodes; ++i) {
    centerline[i][0] = x(2 * i);
    centerline[i][1] = y(2 * i);
    centerline[i][2] = z(2 * i);
    slope[i][0] = x(2 * i + 1);
    slope[i][1] = y(2 * i + 1);
    slope[i][2] = z(2 * i + 1);
  }
}

void
EulerBeamStaticInextensibleADDM::apply_initial_condition_xy()
{
  const real_t h = mesh.get_ds();
  const size_t nodes = mesh.get_nodes();

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

void
EulerBeamStaticInextensibleADDM::apply_initial_condition_xy(EulerBeamMesh& bmesh)
{
  const size_t nodes = bmesh.get_nodes();
  ELFF_ASSERT(nodes == mesh.get_nodes(),
              "Node count of the mesh must match current mesh.\n");
  const auto centerline = bmesh.get_centerline();
  const auto slopes = bmesh.get_slope();

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

void
EulerBeamStaticInextensibleADDM::compute_slopes_collocation()
{
  xp.setZero();
  yp.setZero();
  zp.setZero();
  const real_t h = this->mesh.get_ds();
  const size_t nodes = this->mesh.get_nodes();

  for (size_t ni = 0; ni < nodes; ni++) {
    xp[ni] = x[2 * ni + 1];
    yp[ni] = y[2 * ni + 1];
    zp[ni] = z[2 * ni + 1];
  }

  for (size_t ei = 0; ei < elements; ei++) {
    const std::array<real_t, 4> dH = ELFF::FEM::CubicHermite<real_t>::derivs(0.5, h);

    const size_t edofs[4] = {
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

void
EulerBeamStaticInextensibleADDM::apply_initial_condition_pq()
{
  const size_t nodes = this->mesh.get_nodes();
  compute_slopes_collocation();
  for (size_t ci = 0; ci < nodes + elements; ci++) {
    switch (dimension) {
      case 2: {
        const real_t xprime = xp[ci];
        const real_t yprime = yp[ci];
        const real_t den = std::max(1e-14, sqrt(xprime * xprime + yprime * yprime));
        p[ci] = xprime / den;
        q[ci] = yprime / den;
      } break;
      case 3: {
        const real_t xprime = xp[ci];
        const real_t yprime = yp[ci];
        const real_t zprime = zp[ci];
        const real_t den = std::max(
          1e-14, sqrt(xprime * xprime + yprime * yprime + zprime * zprime));
        p[ci] = xprime / den;
        q[ci] = yprime / den;
        r[ci] = zprime / den;
      } break;
    }
  }
  apply_boundary_condition_pq();
}

void
EulerBeamStaticInextensibleADDM::apply_boundary_condition_pq()
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

void
EulerBeamStaticInextensibleADDM::apply_boundary_condition_lambda()
{
  const size_t nodes = mesh.get_nodes();

  for (size_t bi = 0; bi < 2; ++bi) {
    const EulerBeamBCType bctype = boundary_conditions.type[bi];

    if (bctype != free_bc && bctype != simple_bc) {
      continue;
    }

    size_t ci = 0;
    switch (boundary_conditions.end[bi]) {
      case left:
        ci = 0;
        break;
      case right:
        ci = nodes - 1;
        break;
    }

    lambda_x(ci) = 0.;
    lambda_y(ci) = 0.;
    lambda_z(ci) = 0.;
  }
}

void
EulerBeamStaticInextensibleADDM::update_pq()
{
  const size_t nodes = mesh.get_nodes();
  compute_slopes_collocation();

  for (size_t ci = 0; ci < nodes + elements; ci++) {
    switch (dimension) {
      case 2: {
        const real_t xprime = xp[ci];
        const real_t yprime = yp[ci];
        const real_t X = xprime - lambda_x[ci] / r_penalty;
        const real_t Y = yprime - lambda_y[ci] / r_penalty;
        const real_t norm = std::max(1e-6, sqrt(X * X + Y * Y));
        p[ci] = X / norm;
        q[ci] = Y / norm;
      } break;
      case 3: {
        const real_t xprime = xp[ci];
        const real_t yprime = yp[ci];
        const real_t zprime = zp[ci];
        const real_t X = xprime - lambda_x[ci] / r_penalty;
        const real_t Y = yprime - lambda_y[ci] / r_penalty;
        const real_t Z = zprime - lambda_z[ci] / r_penalty;
        const real_t norm = std::max(1e-6, sqrt(X * X + Y * Y + Z * Z));
        p[ci] = X / norm;
        q[ci] = Y / norm;
        r[ci] = Z / norm;
      } break;
    }
  }

  apply_boundary_condition_pq();
}

void
EulerBeamStaticInextensibleADDM::apply_boundary_condition_A()
{
  A = A_unconstrained;

  const size_t nodes = mesh.get_nodes();

  for (size_t bi = 0; bi < 2; ++bi) {
    const EulerBeamBCEnd bcend = boundary_conditions.end[bi];
    size_t ni = 0;

    std::vector<size_t> idx;
    std::vector<real_t> xvals;
    std::vector<real_t> yvals;
    std::vector<real_t> zvals;

    switch (bcend) {
      case left:
        ni = 0;
        break;
      case right:
        ni = nodes - 1;
        break;
    }

    const EulerBeamBCType bctype = boundary_conditions.type[bi];
    const EulerBeamBCVals bcvals = boundary_conditions.vals[bi];

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
      case point_torque_bc:
        idx = {};
        xvals = {};
        yvals = {};
        zvals = {};
        break;
    }

    for (size_t i = 0; i < xvals.size(); i++) {
      A.row(idx[i]).setZero();
      A.col(idx[i]).setZero();
      A(idx[i], idx[i]) = 1.;
    }
  }
}

void
EulerBeamStaticInextensibleADDM::assemble_A()
{
  const real_t h = this->mesh.get_ds();
  const real_t xi_q[3] = { 0.1127016654, 0.5, 0.8872983346 };
  const real_t w_q[3] = { 0.2777777778, 0.4444444444, 0.2777777778 };

  MatrixXd K4 = MatrixXd::Zero(this->dof, this->dof);
  MatrixXd K2 = MatrixXd::Zero(this->dof, this->dof);

  for (size_t e = 0; e < this->elements; ++e) {
    const size_t edofs[4] = {
      2 * (e + 0) + 0, 2 * (e + 0) + 1, 2 * (e + 1) + 0, 2 * (e + 1) + 1
    };

    real_t K4e[4][4] = { { 0 } };
    real_t K2e[4][4] = { { 0 } };

    for (size_t q = 0; q < 3; ++q) {
      const real_t xi = xi_q[q];
      const real_t w = w_q[q];

      const auto dH = ELFF::FEM::CubicHermite<real_t>::derivs(xi, h);
      const auto ddH = ELFF::FEM::CubicHermite<real_t>::second_derivs(xi, h);

      for (size_t a = 0; a < 4; ++a) {
        for (size_t b = 0; b < 4; ++b) {
          K4e[a][b] += ddH[a] * ddH[b] * w * h;
          K2e[a][b] += dH[a] * dH[b] * w * h;
        }
      }
    }

    for (size_t a = 0; a < 4; ++a) {
      for (size_t b = 0; b < 4; ++b) {
        K4(edofs[a], edofs[b]) += K4e[a][b];
        K2(edofs[a], edofs[b]) += K2e[a][b];
      }
    }
  }

  this->A_unconstrained = this->EI * K4 + this->r_penalty * K2;
  this->A = this->A_unconstrained;
}

void
EulerBeamStaticInextensibleADDM::decompose_A()
{
  llt.compute(A);
}

void
EulerBeamStaticInextensibleADDM::assemble_f(std::array<real_t, 3> load)
{
  f_x.setZero();
  f_y.setZero();
  f_z.setZero();

  const real_t h = this->mesh.get_ds();
  const size_t nodes = this->mesh.get_nodes();

  const real_t xi_q[3] = { 0.1127016654, 0.5, 0.8872983346 };
  const real_t w_q[3] = { 0.2777777778, 0.4444444444, 0.2777777778 };

  for (size_t e = 0; e < elements; ++e) {
    const size_t edofs[4] = {
      2 * (e + 0) + 0, 2 * (e + 0) + 1, 2 * (e + 1) + 0, 2 * (e + 1) + 1
    };
    real_t fxe[4] = { 0, 0, 0, 0 };
    real_t fye[4] = { 0, 0, 0, 0 };
    real_t fze[4] = { 0, 0, 0, 0 };

    for (size_t qi = 0; qi < 3; ++qi) {
      const real_t xi = xi_q[qi];
      const real_t w = w_q[qi];

      const auto L = ELFF::FEM::QuadraticLagrange<real_t>::values(xi);
      const auto H = ELFF::FEM::CubicHermite<real_t>::values(xi, h);
      const auto dH = ELFF::FEM::CubicHermite<real_t>::derivs(xi, h);

      const size_t li = e;
      const size_t mi = nodes + e;
      const size_t ri = e + 1;

      const real_t p_val = L[0] * p[li] + L[1] * p[mi] + L[2] * p[ri];
      const real_t lambda_x_val =
        L[0] * lambda_x[li] + L[1] * lambda_x[mi] + L[2] * lambda_x[ri];
      const real_t q_val = L[0] * q[li] + L[1] * q[mi] + L[2] * q[ri];
      const real_t lambda_y_val =
        L[0] * lambda_y[li] + L[1] * lambda_y[mi] + L[2] * lambda_y[ri];
      const real_t r_val = L[0] * r[li] + L[1] * r[mi] + L[2] * r[ri];
      const real_t lambda_z_val =
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

void
EulerBeamStaticInextensibleADDM::apply_boundary_condition_f()
{
  const size_t nodes = mesh.get_nodes();

  for (size_t bi = 0; bi < 2; ++bi) {
    const EulerBeamBCEnd bcend = boundary_conditions.end[bi];
    size_t ni = 0;

    std::vector<size_t> idx;
    std::vector<real_t> xvals;
    std::vector<real_t> yvals;
    std::vector<real_t> zvals;

    switch (bcend) {
      case left:
        ni = 0;
        break;
      case right:
        ni = nodes - 1;
        break;
    }

    const EulerBeamBCType bctype = boundary_conditions.type[bi];
    const EulerBeamBCVals bcvals = boundary_conditions.vals[bi];

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

    switch (bctype) {
      case point_force_bc:
        for (size_t i = 0; i < xvals.size(); i++) {
          f_x(idx[i]) += xvals[i];
          f_y(idx[i]) += yvals[i];
          f_z(idx[i]) += zvals[i];
        }
        break;
      default:
        for (size_t i = 0; i < xvals.size(); i++) {
          const VectorXd A_col = A_unconstrained.col(idx[i]);

          f_x.noalias() -= A_col * xvals[i];
          f_y.noalias() -= A_col * yvals[i];
          f_z.noalias() -= A_col * zvals[i];

          f_x(idx[i]) = xvals[i];
          f_y(idx[i]) = yvals[i];
          f_z(idx[i]) = zvals[i];
        }
        break;
    }
  }
}

void
EulerBeamStaticInextensibleADDM::update_xy(std::array<real_t, 3> load)
{
  assemble_f(load);
  apply_boundary_condition_f();

  x = llt.solve(f_x);
  y = llt.solve(f_y);
  z = llt.solve(f_z);
}

void
EulerBeamStaticInextensibleADDM::update_multipliers()
{
  const size_t nodes = mesh.get_nodes();
  compute_slopes_collocation();

  for (size_t ci = 0; ci < nodes + elements; ci++) {
    lambda_x[ci] += r_penalty * (p[ci] - xp[ci]);
    lambda_y[ci] += r_penalty * (q[ci] - yp[ci]);
    lambda_z[ci] += r_penalty * (r[ci] - zp[ci]);
  }

  apply_boundary_condition_lambda();
}

bool
EulerBeamStaticInextensibleADDM::is_converged(bool recompute_slopes)
{
  if (recompute_slopes) {
    compute_slopes_collocation();
  }
  const real_t res_p = (p - xp).cwiseAbs().maxCoeff();
  const real_t res_q = (q - yp).cwiseAbs().maxCoeff();
  const real_t res_r = (r - zp).cwiseAbs().maxCoeff();
  const real_t res = std::max({ res_p, res_q, res_r });

  return res < tol_outer;
}

} // namespace Models
} // namespace ELFF
