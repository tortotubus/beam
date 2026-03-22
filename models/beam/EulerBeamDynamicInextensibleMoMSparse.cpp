#include "models/beam/EulerBeamDynamicInextensibleMoMSparse.hpp"

namespace ELFF {
namespace Models {

EulerBeamDynamicInextensibleMoMSparse::EulerBeamDynamicInextensibleMoMSparse(
  real_t length,
  real_t EI,
  real_t mu,
  size_t nodes,
  EulerBeam::EulerBeamBCs bcs,
  real_t r_penalty)
  : EulerBeamStaticInextensibleMoMSparse(length, EI, mu, nodes, bcs, r_penalty)
  , u_prev(VectorXd::Zero(ndof))
  , v_prev(VectorXd::Zero(ndof))
  , a_prev(VectorXd::Zero(ndof))
  , mass(MatrixXd::Zero(ndof, ndof))
{
}

void
EulerBeamDynamicInextensibleMoMSparse::solve(real_t dt,
                                             std::array<real_t, 3> load)
{
  const real_t alpha = 0.0;
  const real_t gamma = 0.5 - alpha;
  const real_t beta = 0.25 * (1 - alpha) * (1 - alpha);
  solve_newmark(dt, load, beta, gamma);
}

void
EulerBeamDynamicInextensibleMoMSparse::solve(
  real_t dt,
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
EulerBeamDynamicInextensibleMoMSparse::solve_newmark(
  real_t dt,
  std::vector<std::array<real_t, 3>> load,
  real_t beta,
  real_t gamma)
{
  ConjugateGradient<SparseMatrix<real_t>, Lower | Upper> solver;

  if (time_iter == 0) {
    VectorXd R0 =
      EulerBeamStaticInextensibleMoMSparse::assemble_residual_template<real_t>(
        u, load);
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
      break;
    } else if (iter_outer == max_iter_outer - 1) {
      ELFF_ABORT(
        "EulerBeamDynamicInextensibleMoMSparse::solve() did not converge.\n");
    }

    solver.setTolerance(tol_inner);
    solver.compute(jacobian);

    VectorXd delta_u = solver.solve(-residual);

    if (solver.info() != Success) {
      ELFF_ABORT(
        "EulerBeamDynamicInextensibleMoMSparse::solve(): linear solve failed.\n");
    }

    u += delta_u;

    S_norm = update_lambda();
  }

  (void) S_norm;

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

  u_prev = u;

  update_mesh();

  time_iter++;
  t += dt;
}

void
EulerBeamDynamicInextensibleMoMSparse::solve_newmark(
  real_t dt,
  std::array<real_t, 3> load,
  real_t beta,
  real_t gamma)
{
  ConjugateGradient<SparseMatrix<real_t>, Lower | Upper> solver;

  if (time_iter == 0) {
    VectorXd R0 =
      EulerBeamStaticInextensibleMoMSparse::assemble_residual_template<real_t>(
        u, load);
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
      break;
    } else if (iter_outer == max_iter_outer - 1) {
      ELFF_ABORT(
        "EulerBeamDynamicInextensibleMoMSparse::solve() did not converge.\n");
    }

    solver.setTolerance(tol_inner);
    solver.compute(jacobian);

    VectorXd delta_u = solver.solve(-residual);

    if (solver.info() != Success) {
      ELFF_ABORT(
        "EulerBeamDynamicInextensibleMoMSparse::solve(): linear solve failed.\n");
    }

    u += delta_u;

    S_norm = update_lambda();
  }

  (void) S_norm;

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

  u_prev = u;

  time_iter++;
  t += dt;
}

void
EulerBeamDynamicInextensibleMoMSparse::apply_initial_condition()
{
  EulerBeamStaticInextensibleMoMSparse::apply_initial_condition();
  u_prev = u;
}

void
EulerBeamDynamicInextensibleMoMSparse::apply_initial_condition(
  EulerBeamMesh& bmesh)
{
  EulerBeamStaticInextensibleMoMSparse::apply_initial_condition(bmesh);
  u_prev = u;
}

void
EulerBeamDynamicInextensibleMoMSparse::assemble_system_newmark(
  real_t dt,
  std::array<real_t, 3> load,
  real_t beta,
  real_t gamma)
{
  static_cast<void>(gamma);

  EulerBeamStaticInextensibleMoMSparse::assemble_system(load);
  add_newmark_inertial_terms(dt, beta);
}

void
EulerBeamDynamicInextensibleMoMSparse::assemble_system_newmark(
  real_t dt,
  std::vector<std::array<real_t, 3>> load,
  real_t beta,
  real_t gamma)
{
  static_cast<void>(gamma);

  EulerBeamStaticInextensibleMoMSparse::assemble_system(load);
  add_newmark_inertial_terms(dt, beta);
}

void
EulerBeamDynamicInextensibleMoMSparse::add_newmark_inertial_terms(real_t dt,
                                                                  real_t beta)
{
  if (!(dt > 0.0)) {
    throw std::runtime_error("Newmark: dt must be > 0");
  }
  if (!(beta > 0.0)) {
    throw std::runtime_error("Newmark: beta must be > 0");
  }

  const real_t inv = 1.0 / (beta * dt * dt);
  const real_t inv_bt = 1.0 / (beta * dt);
  const real_t kappa = (1.0 - 2.0 * beta) / (2.0 * beta);

  for (size_t n = 0; n < nodes; ++n) {
    const size_t ix = offset_x + 2 * n;
    const size_t iy = offset_y + 2 * n;
    const size_t iz = offset_z + 2 * n;
    const real_t w = (n == 0 || n == nodes - 1) ? 0.5 * ds : ds;
    const real_t mass_coeff = mu * w;
    const real_t tangent_coeff = mass_coeff * inv;

    auto add_newmark = [&](size_t i) {
      const real_t a_new =
        inv * (u(i) - u_prev(i)) - inv_bt * v_prev(i) - kappa * a_prev(i);
      residual(i) += mass_coeff * a_new;
      jacobian.coeffRef(i, i) += tangent_coeff;
    };

    add_newmark(ix);
    add_newmark(iy);
    add_newmark(iz);
  }

  jacobian.makeCompressed();
}

void
EulerBeamDynamicInextensibleMoMSparse::update_mesh()
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
