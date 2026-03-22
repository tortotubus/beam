#include "models/beam/EulerBeamStaticInextensibleMoM.hpp"

#include <iostream>

#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/AutoDiff>

namespace ELFF {
namespace Models {

EulerBeamStaticInextensibleMoM::EulerBeamStaticInextensibleMoM(
  real_t length,
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
  , offset_l(ndof_x + ndof_y + ndof_z)
  , r_penalty(r_penalty)
  , max_iter_outer(1000)
  , max_iter_inner(1000)
  , tol_outer(1e-5)
  , tol_inner(1e-5)
  , residual(VectorXd::Zero(ndof))
  , jacobian(MatrixXd::Zero(ndof, ndof))
  , mass()
  , u(VectorXd::Zero(ndof))
  , lambda(VectorXd::Zero(ndof_l))
{
  apply_initial_condition(this->mesh);
}

EulerBeamStaticInextensibleMoM::EulerBeamStaticInextensibleMoM(
  real_t length,
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
  , offset_l(ndof_x + ndof_y + ndof_z)
  , r_penalty(r_penalty)
  , max_iter_outer(1000)
  , max_iter_inner(1000)
  , tol_outer(1e-5)
  , tol_inner(1e-5)
  , residual(VectorXd::Zero(ndof))
  , jacobian(MatrixXd::Zero(ndof, ndof))
  , mass()
  , u(VectorXd::Zero(ndof))
  , lambda(VectorXd::Zero(ndof_l))
{
  apply_initial_condition(this->mesh);
}

EulerBeamStaticInextensibleMoM::~EulerBeamStaticInextensibleMoM() = default;

void
EulerBeamStaticInextensibleMoM::solve()
{
  solve({ 0., 0., 0. });
}

void
EulerBeamStaticInextensibleMoM::solve(std::array<real_t, 3> load)
{
  real_t S_norm = 0;

  LDLT<MatrixXd> solver;

  for (size_t iter_outer = 0; iter_outer < max_iter_outer; iter_outer++) {
    assemble_system(load);
    apply_boundary_conditions();

    real_t res_norm = residual.norm();

    if (res_norm < tol_outer) {
      std::cout << iter_outer << ": ||r|| = " << res_norm;
      std::cout << "\t ||S|| = " << S_norm << std::endl;
      break;
    } else if (iter_outer == max_iter_outer - 1) {
      ELFF_ABORT("EulerBeamStaticInexntensibleMoM::solve() did not converge.\n");
    } else {
      std::cout << iter_outer << ": ||r|| = " << res_norm;
      std::cout << "\t ||S|| = " << S_norm << std::endl;
    }

    solver.compute(jacobian);

    if (solver.info() != Success) {
      ELFF_ABORT("EulerBeamStaticInextensibleMoM::solve(): "
                 "Preconditioner failed.\n");
    }

    VectorXd delta_u = solver.solve(-residual);
    u += delta_u;

    S_norm = update_lambda();
  }

  update_mesh();
}

void
EulerBeamStaticInextensibleMoM::solve(std::vector<std::array<real_t, 3>> load)
{
  real_t S_norm = 0;

  LDLT<MatrixXd> solver;

  for (size_t iter_outer = 0; iter_outer < max_iter_outer; iter_outer++) {
    assemble_system(load);
    apply_boundary_conditions();

    real_t res_norm = residual.norm();

    if (res_norm < tol_outer) {
      std::cout << iter_outer << ": ||r|| = " << res_norm;
      std::cout << "\t ||S|| = " << S_norm << std::endl;
      break;
    } else if (iter_outer == max_iter_outer - 1) {
      ELFF_ABORT("EulerBeamStaticInexntensibleMoM::solve() did not converge.\n");
    } else {
      std::cout << iter_outer << ": ||r|| = " << res_norm;
      std::cout << "\t ||S|| = " << S_norm << std::endl;
    }

    solver.compute(jacobian);

    if (solver.info() != Success) {
      ELFF_ABORT("EulerBeamStaticInextensibleMoM::solve(): "
                 "Preconditioner failed.\n");
    }

    VectorXd delta_u = solver.solve(-residual);
    u += delta_u;

    S_norm = update_lambda();
  }

  update_mesh();
}

void
EulerBeamStaticInextensibleMoM::apply_initial_condition()
{
  for (size_t i = 0; i < nodes; i++) {
    u(offset_x + 2 * i + 0) = ds * i;
    u(offset_x + 2 * i + 1) = 1.;
    u(offset_y + 2 * i + 0) = 0.;
    u(offset_y + 2 * i + 1) = 0.;
    u(offset_z + 2 * i + 0) = 0.;
    u(offset_z + 2 * i + 1) = 0.;
    lambda(i) = 0.;
  }

  update_mesh();
}

void
EulerBeamStaticInextensibleMoM::apply_initial_condition(EulerBeamMesh& bmesh)
{
  if (bmesh.get_nodes() != mesh.get_nodes()) {
    ELFF_ABORT("");
  }

  auto centerline = bmesh.get_centerline();
  auto slope = bmesh.get_slope();

  for (size_t i = 0; i < nodes; i++) {
    u(offset_x + 2 * i + 0) = centerline[i][0];
    u(offset_y + 2 * i + 0) = centerline[i][1];
    u(offset_z + 2 * i + 0) = centerline[i][2];
    u(offset_x + 2 * i + 1) = slope[i][0];
    u(offset_y + 2 * i + 1) = slope[i][1];
    u(offset_z + 2 * i + 1) = slope[i][2];

    lambda(i) = 0.;
  }

  update_mesh();
}

real_t
EulerBeamStaticInextensibleMoM::update_lambda(real_t omega)
{
  // static_cast<void>(omega);

  static constexpr std::array<real_t, 3> xi_q = { 0.1127016654,
                                                  0.5,
                                                  0.8872983346 };
  static constexpr std::array<real_t, 3> w_q = { 0.2777777778,
                                                 0.4444444444,
                                                 0.2777777778 };
  Matrix<real_t, Dynamic, 1> lambda_n = Matrix<real_t, Dynamic, 1>::Zero(ndof_l);

  for (size_t e = 0; e < elements; ++e) {
    size_t elem_nodes[] = { e, e + 1 };
    size_t idx_x[] = { offset_x + 2 * elem_nodes[0],
                       offset_x + 2 * elem_nodes[0] + 1,
                       offset_x + 2 * elem_nodes[1],
                       offset_x + 2 * elem_nodes[1] + 1 };
    size_t idx_y[] = { offset_y + 2 * elem_nodes[0],
                       offset_y + 2 * elem_nodes[0] + 1,
                       offset_y + 2 * elem_nodes[1],
                       offset_y + 2 * elem_nodes[1] + 1 };
    size_t idx_z[] = { offset_z + 2 * elem_nodes[0],
                       offset_z + 2 * elem_nodes[0] + 1,
                       offset_z + 2 * elem_nodes[1],
                       offset_z + 2 * elem_nodes[1] + 1 };
    size_t idx_l[] = { elem_nodes[0], elem_nodes[1] };

    real_t ux[] = { u[idx_x[0]], u[idx_x[1]], u[idx_x[2]], u[idx_x[3]] };
    real_t uy[] = { u[idx_y[0]], u[idx_y[1]], u[idx_y[2]], u[idx_y[3]] };
    real_t uz[] = { u[idx_z[0]], u[idx_z[1]], u[idx_z[2]], u[idx_z[3]] };
    real_t ul[] = { lambda[idx_l[0]], lambda[idx_l[1]] };

    real_t R_loc_l[2] = { 0., 0. };

    for (size_t qi = 0; qi < 3; ++qi) {
      real_t xi = xi_q[qi];
      real_t w = w_q[qi];

      auto dH = CubicHermite<real_t>::derivs(xi, ds);
      auto M = LinearShape<real_t>::values(xi);

      real_t xp = 0;
      real_t yp = 0;
      real_t zp = 0;

      for (size_t i = 0; i < 4; i++) {
        xp += dH[i] * ux[i];
        yp += dH[i] * uy[i];
        zp += dH[i] * uz[i];
      }

      real_t l = 0;
      for (size_t i = 0; i < 2; i++) {
        l += M[i] * ul[i];
      }

      real_t S = xp * xp + yp * yp + zp * zp - 1.0;

      for (size_t a = 0; a < 2; ++a) {
        R_loc_l[a] += S * M[a] * w * ds;
      }
    }

    for (size_t i = 0; i < 2; ++i) {
      lambda_n[idx_l[i]] += R_loc_l[i];
    }
  }

  // lambda = lambda_n;
  lambda = (1.0 - omega) * lambda + omega * lambda_n;


  return lambda_n.norm();
}

void
EulerBeamStaticInextensibleMoM::assemble_residual(std::array<real_t, 3> load)
{
  residual = assemble_residual_template<real_t>(u, load);
}

void
EulerBeamStaticInextensibleMoM::assemble_residual(
  std::vector<std::array<real_t, 3>> load)
{
  residual = assemble_residual_template<real_t>(u, load);
}

void
EulerBeamStaticInextensibleMoM::assemble_system(std::array<real_t, 3> load)
{
  using AD = AutoDiffScalar<VectorXd>;
  using ADVec = Matrix<AD, Dynamic, 1>;

  ADVec x_ad(ndof);
  for (int i = 0; i < int(ndof); ++i) {
    VectorXd seed = VectorXd::Zero(ndof);
    seed(i) = 1.0;
    x_ad(i) = AD(u(i), seed);
  }

  ADVec R_ad = assemble_residual_template<AD>(x_ad, load);

  residual.resize(ndof);
  jacobian.resize(ndof, ndof);
  for (int i = 0; i < int(ndof); ++i) {
    residual(i) = R_ad(i).value();
    jacobian.row(i) = R_ad(i).derivatives().transpose();
  }
}

void
EulerBeamStaticInextensibleMoM::assemble_system(
  std::vector<std::array<real_t, 3>> load)
{
  using AD = AutoDiffScalar<VectorXd>;
  using ADVec = Matrix<AD, Dynamic, 1>;

  ADVec x_ad(ndof);
  for (int i = 0; i < int(ndof); ++i) {
    VectorXd seed = VectorXd::Zero(ndof);
    seed(i) = 1.0;
    x_ad(i) = AD(u(i), seed);
  }

  ADVec R_ad = assemble_residual_template<AD>(x_ad, load);

  residual.resize(ndof);
  jacobian.resize(ndof, ndof);
  for (int i = 0; i < int(ndof); ++i) {
    residual(i) = R_ad(i).value();
    jacobian.row(i) = R_ad(i).derivatives().transpose();
  }
}

void
EulerBeamStaticInextensibleMoM::apply_boundary_conditions()
{
  apply_boundary_conditions(jacobian, residual);
}

void
EulerBeamStaticInextensibleMoM::apply_boundary_conditions(MatrixXd& A,
                                                          VectorXd& R)
{
  (void)A;
  (void)R;

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
        vals = { bcvals.force[0], bcvals.force[1], bcvals.force[2] };
        break;
      case point_torque_bc:
        idx = {
          offset_x + 2 * ni + 1,
          offset_y + 2 * ni + 1,
          offset_z + 2 * ni + 1,
        };
        vals = { bcvals.torque[0], bcvals.torque[1], bcvals.torque[2] };
        break;
      default:
        break;
    }

    switch (bctype) {
      case point_force_bc:
      case point_torque_bc:
        for (size_t i = 0; i < vals.size(); ++i) {
          residual(idx[i]) -= vals[i];
        }
        break;
      default:
        for (size_t i = 0; i < vals.size(); ++i) {
          jacobian.row(idx[i]).setZero();
          jacobian.col(idx[i]).setZero();
          residual(idx[i]) = u(idx[i]) - vals[i];
          jacobian(idx[i], idx[i]) = 1.;
        }
        break;
    }
  }
}

void
EulerBeamStaticInextensibleMoM::update_mesh()
{
  std::vector<std::array<real_t, 3>>& centerline = mesh.get_centerline();
  std::vector<std::array<real_t, 3>>& slope = mesh.get_slope();
  std::vector<real_t>& s = mesh.get_curvilinear_axis();
  (void)s;

  for (size_t i = 0; i < nodes; ++i) {
    centerline[i][0] = u(offset_x + 2 * i + 0);
    centerline[i][1] = u(offset_y + 2 * i + 0);
    centerline[i][2] = u(offset_z + 2 * i + 0);
    slope[i][0] = u(offset_x + 2 * i + 1);
    slope[i][1] = u(offset_y + 2 * i + 1);
    slope[i][2] = u(offset_z + 2 * i + 1);
  }
}

} // namespace Models
} // namespace ELFF
