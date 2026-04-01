#include "elff/models/beam/EulerBeamDynamicInextensibleADDM.hpp"

#include <cmath>

namespace ELFF {
namespace Models {

EulerBeamDynamicInextensibleADDM::EulerBeamDynamicInextensibleADDM(
  real_t length,
  real_t EI,
  real_t mu,
  size_t nodes,
  EulerBeam::EulerBeamBCs bcs,
  real_t r_penalty)
  : EulerBeamStaticInextensibleADDM(length, EI, mu, nodes, bcs, r_penalty)
  , x_prev(VectorXd::Zero(dof))
  , y_prev(VectorXd::Zero(dof))
  , z_prev(VectorXd::Zero(dof))
  , vx_prev(VectorXd::Zero(dof))
  , vy_prev(VectorXd::Zero(dof))
  , vz_prev(VectorXd::Zero(dof))
  , ax_prev(VectorXd::Zero(dof))
  , ay_prev(VectorXd::Zero(dof))
  , az_prev(VectorXd::Zero(dof))
  , mass(MatrixXd::Zero(dof, dof))
  , A_static_unconstrained(A_unconstrained)
  , load_prev({ 0., 0., 0. })
{
  max_outer = 1000;
  tol_outer = 1e-5;
  x_prev = x;
  y_prev = y;
  z_prev = z;
  assemble_mass_matrix();
}

void
EulerBeamDynamicInextensibleADDM::solve(real_t dt, std::array<real_t, 3> load)
{
  const real_t alpha = 0.0;
  const real_t gamma = 0.5 - alpha;
  const real_t beta = 0.25 * (1.0 - alpha) * (1.0 - alpha);
  solve_newmark(dt, load, beta, gamma);
}

void
EulerBeamDynamicInextensibleADDM::solve(real_t dt,
                                        std::vector<std::array<real_t, 3>> load)
{
  ELFF_ASSERT(load.size() == mesh.get_nodes(),
              "Size of load vector must equal number of nodes.");

  const real_t alpha = 0.0;
  const real_t gamma = 0.5 - alpha;
  const real_t beta = 0.25 * (1.0 - alpha) * (1.0 - alpha);
  solve_newmark(dt, load, beta, gamma);
}

void
EulerBeamDynamicInextensibleADDM::solve_newmark(real_t dt,
                                                std::array<real_t, 3> load,
                                                real_t beta,
                                                real_t gamma)
{
  if (!(dt > 0.0)) {
    throw std::runtime_error("Newmark: dt must be > 0");
  }
  if (!(beta > 0.0)) {
    throw std::runtime_error("Newmark: beta must be > 0");
  }
  if (!(gamma > 0.0)) {
    throw std::runtime_error("Newmark: gamma must be > 0");
  }

  x = x_prev;
  y = y_prev;
  z = z_prev;
  lambda_x.setZero();
  lambda_y.setZero();
  lambda_z.setZero();
  apply_initial_condition_pq();
  apply_boundary_condition_lambda();

  const VectorXd x_old = x_prev;
  const VectorXd y_old = y_prev;
  const VectorXd z_old = z_prev;
  real_t final_rel_update = 0.0;
  real_t final_max_pq_error = 0.0;
  real_t final_max_xy_error = 0.0;
  bool converged = false;
  size_t iter;
  for (iter = 0; iter < max_outer; ++iter) {
    const VectorXd x_iter_prev = x;
    const VectorXd y_iter_prev = y;
    update_pq();
    assemble_system_newmark(load, dt, beta, gamma);
    x = llt.solve(f_x);
    y = llt.solve(f_y);
    z = llt.solve(f_z);
    update_multipliers();
    final_rel_update =
      compute_relative_xy_update(x_iter_prev, x, y_iter_prev, y);
    final_max_pq_error = compute_max_pq_error();
    final_max_xy_error =
      compute_max_xy_update(x_iter_prev, x, y_iter_prev, y);

    if (final_rel_update < tol_outer) {
      converged = true;
      break;
    }
  }

  if (!converged) {
    ELFF_WARNING("EulerBeamDynamicInextensibleADDM::solve_newmark() final "
                 "relative xy update = "
                 << final_rel_update << " at step " << time_iter << " after "
                 << iter << " iterations");
    // ELFF_ABORT("EulerBeamDynamicInextensibleADDM::solve_newmark() did not "
    //            "converge.\n");
  }

  ELFF_LOG(time_iter << "\t" << final_max_pq_error << "\t"
                     << final_max_xy_error << "\t" << iter);

  update_newmark_state_component(x_old, x, vx_prev, ax_prev, dt, beta, gamma);
  update_newmark_state_component(y_old, y, vy_prev, ay_prev, dt, beta, gamma);
  update_newmark_state_component(z_old, z, vz_prev, az_prev, dt, beta, gamma);
  apply_dynamic_state_boundary_conditions();

  x_prev = x;
  y_prev = y;
  z_prev = z;
  load_prev = load;

  update_mesh();

  ++time_iter;
  t += dt;
}

void
EulerBeamDynamicInextensibleADDM::solve_newmark(
  real_t dt,
  std::vector<std::array<real_t, 3>> load,
  real_t beta,
  real_t gamma)
{
  ELFF_ASSERT(load.size() == mesh.get_nodes(),
              "Size of load vector must equal number of nodes.");

  if (!(dt > 0.0)) {
    throw std::runtime_error("Newmark: dt must be > 0");
  }
  if (!(beta > 0.0)) {
    throw std::runtime_error("Newmark: beta must be > 0");
  }
  if (!(gamma > 0.0)) {
    throw std::runtime_error("Newmark: gamma must be > 0");
  }

  x = x_prev;
  y = y_prev;
  z = z_prev;
  lambda_x.setZero();
  lambda_y.setZero();
  lambda_z.setZero();
  apply_initial_condition_pq();
  apply_boundary_condition_lambda();

  const VectorXd x_old = x_prev;
  const VectorXd y_old = y_prev;
  const VectorXd z_old = z_prev;
  real_t final_rel_update = 0.0;
  real_t final_max_pq_error = 0.0;
  real_t final_max_xy_error = 0.0;
  bool converged = false;
  size_t iter;

  for (iter = 0; iter < max_outer; ++iter) {
    const VectorXd x_iter_prev = x;
    const VectorXd y_iter_prev = y;
    update_pq();
    assemble_system_newmark(load, dt, beta, gamma);
    x = llt.solve(f_x);
    y = llt.solve(f_y);
    z = llt.solve(f_z);
    update_multipliers();
    final_rel_update =
      compute_relative_xy_update(x_iter_prev, x, y_iter_prev, y);
    final_max_pq_error = compute_max_pq_error();
    final_max_xy_error =
      compute_max_xy_update(x_iter_prev, x, y_iter_prev, y);

    if (final_rel_update < tol_outer) {
      converged = true;
      break;
    }
  }

  if (!converged) {
    ELFF_WARNING("EulerBeamDynamicInextensibleADDM::solve_newmark() final "
                 "relative xy update = "
                 << final_rel_update << " at step " << time_iter << " after "
                 << iter << " iterations");
    ELFF_ABORT("EulerBeamDynamicInextensibleADDM::solve_newmark() did not "
               "converge.\n");
  }

  ELFF_LOG(time_iter << "\t" << final_max_pq_error << "\t"
                     << final_max_xy_error << "\t" << iter);

  update_newmark_state_component(x_old, x, vx_prev, ax_prev, dt, beta, gamma);
  update_newmark_state_component(y_old, y, vy_prev, ay_prev, dt, beta, gamma);
  update_newmark_state_component(z_old, z, vz_prev, az_prev, dt, beta, gamma);
  apply_dynamic_state_boundary_conditions();

  x_prev = x;
  y_prev = y;
  z_prev = z;
  load_prev = { 0., 0., 0. };

  update_mesh();

  ++time_iter;
  t += dt;
}

void
EulerBeamDynamicInextensibleADDM::apply_initial_condition()
{
  EulerBeamStaticInextensibleADDM::apply_initial_condition();
  x_prev = x;
  y_prev = y;
  z_prev = z;
  vx_prev.setZero();
  vy_prev.setZero();
  vz_prev.setZero();
  ax_prev.setZero();
  ay_prev.setZero();
  az_prev.setZero();
  load_prev = { 0., 0., 0. };
  update_mesh();
}

void
EulerBeamDynamicInextensibleADDM::apply_initial_condition(EulerBeamMesh& bmesh)
{
  EulerBeamStaticInextensibleADDM::apply_initial_condition(bmesh);
  const auto velocities = bmesh.get_centerline_velocity();
  x_prev = x;
  y_prev = y;
  z_prev = z;
  vx_prev.setZero();
  vy_prev.setZero();
  vz_prev.setZero();
  for (size_t i = 0; i < mesh.get_nodes(); ++i) {
    vx_prev(2 * i + 0) = velocities[i][0];
    vy_prev(2 * i + 0) = velocities[i][1];
    vz_prev(2 * i + 0) = velocities[i][2];
  }
  ax_prev.setZero();
  ay_prev.setZero();
  az_prev.setZero();
  load_prev = { 0., 0., 0. };
  update_mesh();
}

void
EulerBeamDynamicInextensibleADDM::assemble_mass_matrix()
{
  mass.setZero();

  const real_t h = mesh.get_ds();
  const size_t nodes = mesh.get_nodes();

  for (size_t ni = 0; ni < nodes; ++ni) {
    const real_t w = (ni == 0 || ni == nodes - 1) ? 0.5 * h : h;
    mass(2 * ni + 0, 2 * ni + 0) = mu * w;
  }
}

void
EulerBeamDynamicInextensibleADDM::assemble_system_newmark(
  std::array<real_t, 3> load,
  real_t dt,
  real_t beta,
  real_t gamma)
{
  static_cast<void>(gamma);

  EulerBeamStaticInextensibleADDM::assemble_f(load);

  const real_t coeff = 1.0 / (beta * dt * dt);
  const real_t inv_bt = 1.0 / (beta * dt);
  const real_t kappa = (1.0 - 2.0 * beta) / (2.0 * beta);

  const VectorXd x_predict =
    coeff * x_prev + inv_bt * vx_prev + kappa * ax_prev;
  const VectorXd y_predict =
    coeff * y_prev + inv_bt * vy_prev + kappa * ay_prev;
  const VectorXd z_predict =
    coeff * z_prev + inv_bt * vz_prev + kappa * az_prev;

  f_x.noalias() += mass * x_predict;
  f_y.noalias() += mass * y_predict;
  f_z.noalias() += mass * z_predict;

  A_unconstrained = A_static_unconstrained + coeff * mass;
  A = A_unconstrained;

  apply_boundary_condition_A();
  apply_boundary_condition_f();
  decompose_A();
}

void
EulerBeamDynamicInextensibleADDM::assemble_system_newmark(
  const std::vector<std::array<real_t, 3>>& load,
  real_t dt,
  real_t beta,
  real_t gamma)
{
  static_cast<void>(gamma);

  assemble_f_nodal(load);

  const real_t coeff = 1.0 / (beta * dt * dt);
  const real_t inv_bt = 1.0 / (beta * dt);
  const real_t kappa = (1.0 - 2.0 * beta) / (2.0 * beta);

  const VectorXd x_predict =
    coeff * x_prev + inv_bt * vx_prev + kappa * ax_prev;
  const VectorXd y_predict =
    coeff * y_prev + inv_bt * vy_prev + kappa * ay_prev;
  const VectorXd z_predict =
    coeff * z_prev + inv_bt * vz_prev + kappa * az_prev;

  f_x.noalias() += mass * x_predict;
  f_y.noalias() += mass * y_predict;
  f_z.noalias() += mass * z_predict;

  A_unconstrained = A_static_unconstrained + coeff * mass;
  A = A_unconstrained;

  apply_boundary_condition_A();
  apply_boundary_condition_f();
  decompose_A();
}

void
EulerBeamDynamicInextensibleADDM::assemble_f_nodal(
  const std::vector<std::array<real_t, 3>>& load)
{
  ELFF_ASSERT(load.size() == mesh.get_nodes(),
              "Size of load vector must equal number of nodes.");

  f_x.setZero();
  f_y.setZero();
  f_z.setZero();

  const real_t h = mesh.get_ds();
  const size_t nodes = mesh.get_nodes();

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
      const auto M = ELFF::FEM::LinearShape<real_t>::values(xi);
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

      const real_t fx_q = M[0] * load[e][0] + M[1] * load[e + 1][0];
      const real_t fy_q = M[0] * load[e][1] + M[1] * load[e + 1][1];
      const real_t fz_q = M[0] * load[e][2] + M[1] * load[e + 1][2];

      for (size_t a = 0; a < 4; ++a) {
        fxe[a] += (lambda_x_val + r_penalty * p_val) * dH[a] * w * h;
        fye[a] += (lambda_y_val + r_penalty * q_val) * dH[a] * w * h;
        fze[a] += (lambda_z_val + r_penalty * r_val) * dH[a] * w * h;
        fxe[a] += fx_q * H[a] * w * h;
        fye[a] += fy_q * H[a] * w * h;
        fze[a] += fz_q * H[a] * w * h;
      }
    }

    for (size_t a = 0; a < 4; ++a) {
      f_x(edofs[a]) += fxe[a];
      f_y(edofs[a]) += fye[a];
      f_z(edofs[a]) += fze[a];
    }
  }

  for (size_t bi = 0; bi < 2; ++bi) {
    const EulerBeamBCEnd bcend = boundary_conditions.end[bi];
    const EulerBeamBCType bctype = boundary_conditions.type[bi];
    const EulerBeamBCVals bcvals = boundary_conditions.vals[bi];

    size_t ni = 0;
    switch (bcend) {
      case left:
        ni = 0;
        break;
      case right:
        ni = nodes - 1;
        break;
    }

    if (bctype == point_force_bc) {
      f_x(2 * ni + 0) += bcvals.force[0];
      f_y(2 * ni + 0) += bcvals.force[1];
      f_z(2 * ni + 0) += bcvals.force[2];
    } else if (bctype == point_torque_bc) {
      f_x(2 * ni + 1) += bcvals.torque[0];
      f_y(2 * ni + 1) += bcvals.torque[1];
      f_z(2 * ni + 1) += bcvals.torque[2];
    }
  }
}

void
EulerBeamDynamicInextensibleADDM::initialize_newmark_acceleration(
  std::array<real_t, 3> load)
{
  update_pq();
  EulerBeamStaticInextensibleADDM::assemble_f(load);

  VectorXd rhs_x = f_x - A_static_unconstrained * x;
  VectorXd rhs_y = f_y - A_static_unconstrained * y;
  VectorXd rhs_z = f_z - A_static_unconstrained * z;

  std::vector<size_t> idx;
  std::vector<real_t> xvals, yvals, zvals;
  collect_boundary_dofs(idx, xvals, yvals, zvals);

  MatrixXd mass_bc = mass;
  for (const size_t d : idx) {
    mass_bc.row(d).setZero();
    mass_bc.col(d).setZero();
    mass_bc(d, d) = 1.0;
    rhs_x(d) = 0.0;
    rhs_y(d) = 0.0;
    rhs_z(d) = 0.0;
  }

  LLT<MatrixXd> mass_solver;
  mass_solver.compute(mass_bc);
  if (mass_solver.info() != Success) {
    ELFF_ABORT("EulerBeamDynamicInextensibleADDM::initialize_newmark_"
               "acceleration(): mass factorization failed.\n");
  }

  ax_prev = mass_solver.solve(rhs_x);
  ay_prev = mass_solver.solve(rhs_y);
  az_prev = mass_solver.solve(rhs_z);
  apply_dynamic_state_boundary_conditions();
}

void
EulerBeamDynamicInextensibleADDM::initialize_newmark_acceleration(
  const std::vector<std::array<real_t, 3>>& load)
{
  update_pq();
  assemble_f_nodal(load);

  VectorXd rhs_x = f_x - A_static_unconstrained * x;
  VectorXd rhs_y = f_y - A_static_unconstrained * y;
  VectorXd rhs_z = f_z - A_static_unconstrained * z;

  std::vector<size_t> idx;
  std::vector<real_t> xvals, yvals, zvals;
  collect_boundary_dofs(idx, xvals, yvals, zvals);

  MatrixXd mass_bc = mass;
  for (const size_t d : idx) {
    mass_bc.row(d).setZero();
    mass_bc.col(d).setZero();
    mass_bc(d, d) = 1.0;
    rhs_x(d) = 0.0;
    rhs_y(d) = 0.0;
    rhs_z(d) = 0.0;
  }

  LLT<MatrixXd> mass_solver;
  mass_solver.compute(mass_bc);
  if (mass_solver.info() != Success) {
    ELFF_ABORT("EulerBeamDynamicInextensibleADDM::initialize_newmark_"
               "acceleration(): mass factorization failed.\n");
  }

  ax_prev = mass_solver.solve(rhs_x);
  ay_prev = mass_solver.solve(rhs_y);
  az_prev = mass_solver.solve(rhs_z);
  apply_dynamic_state_boundary_conditions();
}

void
EulerBeamDynamicInextensibleADDM::apply_dynamic_state_boundary_conditions()
{
  const size_t nodes = mesh.get_nodes();

  for (size_t bi = 0; bi < 2; ++bi) {
    const EulerBeamBCType bctype = boundary_conditions.type[bi];

    if (bctype != simple_bc && bctype != clamped_bc) {
      continue;
    }

    size_t ni = 0;
    switch (boundary_conditions.end[bi]) {
      case left:
        ni = 0;
        break;
      case right:
        ni = nodes - 1;
        break;
    }

    vx_prev(2 * ni + 0) = 0.;
    vy_prev(2 * ni + 0) = 0.;
    vz_prev(2 * ni + 0) = 0.;
    ax_prev(2 * ni + 0) = 0.;
    ay_prev(2 * ni + 0) = 0.;
    az_prev(2 * ni + 0) = 0.;

    if (bctype == clamped_bc) {
      vx_prev(2 * ni + 1) = 0.;
      vy_prev(2 * ni + 1) = 0.;
      vz_prev(2 * ni + 1) = 0.;
      ax_prev(2 * ni + 1) = 0.;
      ay_prev(2 * ni + 1) = 0.;
      az_prev(2 * ni + 1) = 0.;
    }
  }
}

real_t
EulerBeamDynamicInextensibleADDM::compute_relative_xy_update(
  const VectorXd& x_old_iter,
  const VectorXd& x_new_iter,
  const VectorXd& y_old_iter,
  const VectorXd& y_new_iter) const
{
  const VectorXd dx = x_new_iter - x_old_iter;
  const VectorXd dy = y_new_iter - y_old_iter;

  const real_t numerator = std::sqrt(dx.squaredNorm() + dy.squaredNorm());
  const real_t denominator = std::max<real_t>(
    1e-14, std::sqrt(x_new_iter.squaredNorm() + y_new_iter.squaredNorm()));

  return numerator / denominator;
}

real_t
EulerBeamDynamicInextensibleADDM::compute_max_pq_error() const
{
  const real_t res_p = (p - xp).cwiseAbs().maxCoeff();
  const real_t res_q = (q - yp).cwiseAbs().maxCoeff();
  const real_t res_r = (r - zp).cwiseAbs().maxCoeff();
  return std::max(res_p, std::max(res_q, res_r));
}

real_t
EulerBeamDynamicInextensibleADDM::compute_max_xy_update(
  const VectorXd& x_old_iter,
  const VectorXd& x_new_iter,
  const VectorXd& y_old_iter,
  const VectorXd& y_new_iter) const
{
  const real_t max_dx = (x_new_iter - x_old_iter).cwiseAbs().maxCoeff();
  const real_t max_dy = (y_new_iter - y_old_iter).cwiseAbs().maxCoeff();
  return std::max(max_dx, max_dy);
}

void
EulerBeamDynamicInextensibleADDM::update_newmark_state_component(
  const VectorXd& u_old,
  const VectorXd& u_new,
  VectorXd& v_hist,
  VectorXd& a_hist,
  real_t dt,
  real_t beta,
  real_t gamma)
{
  const VectorXd v_old = v_hist;
  const VectorXd a_old = a_hist;
  const real_t inv = 1.0 / (beta * dt * dt);
  const real_t kappa = (1.0 - 2.0 * beta) / (2.0 * beta);

  a_hist = inv * (u_new - u_old - dt * v_old) - kappa * a_old;
  v_hist = v_old + dt * ((1.0 - gamma) * a_old + gamma * a_hist);
}

void
EulerBeamDynamicInextensibleADDM::update_mesh()
{
  EulerBeamStaticInextensibleADDM::update_mesh();

  std::vector<std::array<real_t, 3>>& velocity = mesh.get_centerline_velocity();
  for (size_t i = 0; i < mesh.get_nodes(); ++i) {
    velocity[i][0] = vx_prev(2 * i + 0);
    velocity[i][1] = vy_prev(2 * i + 0);
    velocity[i][2] = vz_prev(2 * i + 0);
  }
}

} // namespace Models
} // namespace ELFF
