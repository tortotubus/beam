#include "models/beam/EulerBeamDynamicInextensibleMoM.hpp"

namespace ELFF {
namespace Models {

EulerBeamDynamicInextensibleMoM::EulerBeamDynamicInextensibleMoM(
  real_t length,
  real_t EI,
  real_t mu,
  size_t nodes,
  EulerBeam::EulerBeamBCs bcs,
  real_t r_penalty)
  : EulerBeamStaticInextensibleMoM(length, EI, mu, nodes, bcs, r_penalty)
  , u_prev(VectorXd::Zero(ndof))
  , v_prev(VectorXd::Zero(ndof))
  , a_prev(VectorXd::Zero(ndof))
  , mass(MatrixXd::Zero(ndof, ndof))
{
}

void
EulerBeamDynamicInextensibleMoM::solve(real_t dt,
                                       std::array<real_t, 3> load)
{
  const real_t alpha = 0.0;
  const real_t gamma = 0.5 - alpha;
  const real_t beta = 0.25 * (1 - alpha) * (1 - alpha);
  solve_newmark(dt, load, beta, gamma);
}

void
EulerBeamDynamicInextensibleMoM::solve(real_t dt,
                                       std::vector<std::array<real_t, 3>> load)
{
  ELFF_ASSERT(load.size() == nodes,
              "Size of load vector must equal number of nodes.");
  const real_t alpha = 0.0;
  const real_t gamma = 0.5 - alpha;
  const real_t beta = 0.25 * (1 - alpha) * (1 - alpha);
  solve_newmark(dt, load, beta, gamma);
}

void
EulerBeamDynamicInextensibleMoM::solve_newmark(
  real_t dt,
  std::vector<std::array<real_t, 3>> load,
  real_t beta,
  real_t gamma)
{

  // LDLT<MatrixXd> solver;
  // solver.setTolerance(tol_inner);
  ConjugateGradient<MatrixXd, Upper | Lower> solver;

  if (time_iter == 0) {
    VectorXd R0 =
      EulerBeamStaticInextensibleMoM ::assemble_residual_template<real_t>(u,
                                                                          load);
    apply_boundary_conditions();
    for (size_t n = 0; n < nodes; ++n) {
      size_t ix = offset_x + 2 * n;
      size_t iy = offset_y + 2 * n;
      size_t iz = offset_z + 2 * n;
      const real_t w = (n == 0 || n == nodes - 1) ? 0.5 * ds : ds;
      a_prev(ix) = (-R0(ix)) / (mu * w);
      a_prev(iy) = (-R0(iy)) / (mu * w);
      a_prev(iz) = (-R0(iz)) / (mu * w);
    }
  }

  u_prev = u;

  real_t S_norm = 0;

  for (size_t iter_outer = 0; iter_outer < max_iter_outer; iter_outer++) {
    assemble_system_newmark(dt, load, beta, gamma);
    apply_boundary_conditions();

    real_t res_norm = residual.norm();

    if (res_norm < tol_outer) {
      // std::cout << time_iter << " : " << iter_outer << " : ";
      // std::cout << "||r|| = " << res_norm << "\t";
      // std::cout << "||S|| = " << S_norm << std::endl;
      break;
    } else if (iter_outer == max_iter_outer - 1) {
      ELFF_ABORT(
        "EulerBeamDynamicInextensibleMoM::solve() did not converge.\n");
    } else {
      // std::cout << time_iter << " : " << iter_outer << " : ";
      // std::cout << "||r|| = " << res_norm << "\t";
      // std::cout << "||S|| = " << S_norm << std::endl;
    }

    solver.compute(jacobian);
    VectorXd delta_u = solver.solve(-residual);
    u += delta_u;

    S_norm = update_lambda();
  }

  size_t nodes = mesh.get_nodes();

  for (size_t ni = 0; ni < nodes; ++ni) {
    size_t ix = offset_x + 2 * ni;
    size_t iy = offset_y + 2 * ni;
    size_t iz = offset_z + 2 * ni;

    auto upd = [&](size_t i) {
      real_t a_new = (u(i) - u_prev(i) - dt * v_prev(i)) / (beta * dt * dt) -
                     ((1.0 - 2.0 * beta) / (2.0 * beta)) * a_prev(i);
      real_t v_new =
        v_prev(i) + dt * ((1.0 - gamma) * a_prev(i) + gamma * a_new);
      a_prev(i) = a_new;
      v_prev(i) = v_new;
    };

    upd(ix);
    upd(iy);
    upd(iz);
  }

  u_prev = u;           // store new n

  update_mesh();

  time_iter++;
  t += dt;
}

void
EulerBeamDynamicInextensibleMoM::solve_newmark(real_t dt,
                                               std::array<real_t, 3> load,
                                               real_t beta,
                                               real_t gamma)
{

  // LDLT<MatrixXd> solver;
  // solver.setTolerance(tol_inner);
  ConjugateGradient<MatrixXd, Upper | Lower> solver;

  if (time_iter == 0) {
    VectorXd R0 =
      EulerBeamStaticInextensibleMoM ::assemble_residual_template<real_t>(u,
                                                                          load);
    apply_boundary_conditions();
    for (size_t n = 0; n < nodes; ++n) {
      size_t ix = offset_x + 2 * n;
      size_t iy = offset_y + 2 * n;
      size_t iz = offset_z + 2 * n;
      const real_t w = (n == 0 || n == nodes - 1) ? 0.5 * ds : ds;
      a_prev(ix) = (-R0(ix)) / (mu * w);
      a_prev(iy) = (-R0(iy)) / (mu * w);
      a_prev(iz) = (-R0(iz)) / (mu * w);
    }
  }

  u_prev = u;

  real_t S_norm = 0;

  for (size_t iter_outer = 0; iter_outer < max_iter_outer; iter_outer++) {
    assemble_system_newmark(dt, load, beta, gamma);
    apply_boundary_conditions();

    real_t res_norm = residual.norm();

    if (res_norm < tol_outer) {
      // std::cout << time_iter << " : " << iter_outer << " : ";
      // std::cout << "||r|| = " << res_norm << "\t";
      // std::cout << "||S|| = " << S_norm << std::endl;
      break;
    } else if (iter_outer == max_iter_outer - 1) {
      ELFF_ABORT(
        "EulerBeamDynamicInextensibleMoM::solve() did not converge.\n");
    } else {
      // std::cout << time_iter << " : " << iter_outer << " : ";
      // std::cout << "||r|| = " << res_norm << "\t";
      // std::cout << "||S|| = " << S_norm << std::endl;
    }

    solver.compute(jacobian);
    VectorXd delta_u = solver.solve(-residual);
    u += delta_u;

    S_norm = update_lambda();
  }

  update_mesh();

  size_t nodes = mesh.get_nodes();

  for (size_t ni = 0; ni < nodes; ++ni) {
    size_t ix = offset_x + 2 * ni;
    size_t iy = offset_y + 2 * ni;
    size_t iz = offset_z + 2 * ni;

    auto upd = [&](size_t i) {
      real_t a_new = (u(i) - u_prev(i) - dt * v_prev(i)) / (beta * dt * dt) -
                     ((1.0 - 2.0 * beta) / (2.0 * beta)) * a_prev(i);
      real_t v_new =
        v_prev(i) + dt * ((1.0 - gamma) * a_prev(i) + gamma * a_new);
      a_prev(i) = a_new;
      v_prev(i) = v_new;
    };

    upd(ix);
    upd(iy);
    upd(iz);
  }

  u_prev = u;           // store new n

  time_iter++;
  t += dt;
}

void
EulerBeamDynamicInextensibleMoM::apply_initial_condition()
{
  EulerBeamStaticInextensibleMoM::apply_initial_condition();
  u_prev = u;
}

void
EulerBeamDynamicInextensibleMoM::apply_initial_condition(EulerBeamMesh& bmesh)
{
  EulerBeamStaticInextensibleMoM::apply_initial_condition(bmesh);
  u_prev = u;
}

void
EulerBeamDynamicInextensibleMoM::assemble_system_newmark(
  real_t dt,
  std::array<real_t, 3> load,
  real_t beta,
  real_t gamma)
{
  using AD = AutoDiffScalar<VectorXd>;
  using ADVec = Matrix<AD, Dynamic, 1>;

  ADVec x_ad(ndof);

  for (int i = 0; i < int(ndof); ++i) {
    VectorXd seed = VectorXd::Zero(ndof);
    seed(i) = 1.0;
    x_ad(i) = AD(u(i), seed);
  }

  ADVec R_ad = assemble_residual_newmark<AD>(x_ad, dt, load, beta, gamma);

  residual.resize(ndof);
  jacobian.resize(ndof, ndof);
  for (int i = 0; i < int(ndof); ++i) {
    residual(i) = R_ad(i).value();
    jacobian.row(i) = R_ad(i).derivatives().transpose();
  }
}

void
EulerBeamDynamicInextensibleMoM::assemble_system_newmark(
  real_t dt,
  std::vector<std::array<real_t, 3>> load,
  real_t beta,
  real_t gamma)
{
  using AD = AutoDiffScalar<VectorXd>;
  using ADVec = Matrix<AD, Dynamic, 1>;

  ADVec x_ad(ndof);

  for (int i = 0; i < int(ndof); ++i) {
    VectorXd seed = VectorXd::Zero(ndof);
    seed(i) = 1.0;
    x_ad(i) = AD(u(i), seed);
  }

  ADVec R_ad = assemble_residual_newmark<AD>(x_ad, dt, load, beta, gamma);

  residual.resize(ndof);
  jacobian.resize(ndof, ndof);
  for (int i = 0; i < int(ndof); ++i) {
    residual(i) = R_ad(i).value();
    jacobian.row(i) = R_ad(i).derivatives().transpose();
  }
}

void
EulerBeamDynamicInextensibleMoM::update_mesh()
{

  std::vector<std::array<real_t, 3>>& centerline = mesh.get_centerline();
  std::vector<std::array<real_t, 3>>& velocity = mesh.get_centerline_velocity();
  std::vector<std::array<real_t, 3>>& slope = mesh.get_slope();
  std::vector<real_t>& s = mesh.get_curvilinear_axis();

  for (size_t i = 0; i < nodes; ++i) {
    centerline[i][0] = u(offset_x + 2 * i + 0);
    centerline[i][1] = u(offset_y + 2 * i + 0);
    centerline[i][2] = u(offset_z + 2 * i + 0);
    slope[i][0] = u(offset_x + 2 * i + 1);
    slope[i][1] = u(offset_y + 2 * i + 1);
    slope[i][2] = u(offset_z + 2 * i + 1);
    velocity[i][0] = v_prev(offset_x + 2 * i + 0);
    velocity[i][1] = v_prev(offset_y + 2 * i + 0);
    velocity[i][2] = v_prev(offset_z + 2 * i + 0);
  }
}

} // namespace Models
} // namespace ELFF
