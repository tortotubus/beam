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
  , load_prev({ 0., 0., 0. })
  , nodal_load_prev()
  , have_prev_uniform_load(false)
  , have_prev_nodal_load(false)
{
  max_outer = 2000;
  tol_outer = 1e-7;
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
  ELFF_ASSERT(std::abs(beta - 0.25) < 1e-12 && std::abs(gamma - 0.5) < 1e-12,
              "EulerBeamDynamicInextensibleADDM implements the average-"
              "acceleration Newmark variant used in the paper.\n");

  if (!have_prev_uniform_load) {
    load_prev = load;
    have_prev_uniform_load = true;
  }
  have_prev_nodal_load = false;

  x = x_prev;
  y = y_prev;
  z = z_prev;

  const VectorXd x_old = x_prev;
  const VectorXd y_old = y_prev;
  const VectorXd z_old = z_prev;
  real_t final_rel_update = 0.0;
  real_t final_max_pq_error = 0.0;
  real_t final_max_state_error = 0.0;
  bool converged = false;
  size_t iter;
  for (iter = 0; iter < max_outer; ++iter) {
    const VectorXd x_iter_prev = x;
    const VectorXd y_iter_prev = y;
    const VectorXd z_iter_prev = z;
    update_pq();
    assemble_system_newmark(load, dt);
    x = llt.solve(f_x);
    y = llt.solve(f_y);
    z = llt.solve(f_z);
    update_multipliers();
    final_rel_update = compute_relative_state_update(
      x_iter_prev, x, y_iter_prev, y, z_iter_prev, z);
    final_max_pq_error = compute_max_pq_error();
    final_max_state_error = compute_max_state_update(
      x_iter_prev, x, y_iter_prev, y, z_iter_prev, z);

    if (final_rel_update < tol_outer) {
      converged = true;
      break;
    }
  }

  if (!converged) {
    ELFF_WARNING("EulerBeamDynamicInextensibleADDM::solve_newmark() final "
                 "relative state update = "
                 << final_rel_update << " at step " << time_iter << " after "
                 << iter << " iterations");
    // ELFF_ABORT("EulerBeamDynamicInextensibleADDM::solve_newmark() did not "
    //            "converge.\n");
  }

  ELFF_LOG(time_iter << "\t" << final_rel_update << "\t" << iter);

  update_average_acceleration_state_component(x_old, x, vx_prev, ax_prev, dt);
  update_average_acceleration_state_component(y_old, y, vy_prev, ay_prev, dt);
  update_average_acceleration_state_component(z_old, z, vz_prev, az_prev, dt);
  apply_dynamic_state_boundary_conditions();

  x_prev = x;
  y_prev = y;
  z_prev = z;
  load_prev = load;
  have_prev_uniform_load = true;

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
  ELFF_ASSERT(std::abs(beta - 0.25) < 1e-12 && std::abs(gamma - 0.5) < 1e-12,
              "EulerBeamDynamicInextensibleADDM implements the average-"
              "acceleration Newmark variant used in the paper.\n");

  if (!have_prev_nodal_load) {
    nodal_load_prev = load;
    have_prev_nodal_load = true;
  }
  have_prev_uniform_load = false;

  x = x_prev;
  y = y_prev;
  z = z_prev;

  const VectorXd x_old = x_prev;
  const VectorXd y_old = y_prev;
  const VectorXd z_old = z_prev;
  real_t final_rel_update = 0.0;
  real_t final_max_pq_error = 0.0;
  real_t final_max_state_error = 0.0;
  bool converged = false;
  size_t iter;

  for (iter = 0; iter < max_outer; ++iter) {
    const VectorXd x_iter_prev = x;
    const VectorXd y_iter_prev = y;
    const VectorXd z_iter_prev = z;
    update_pq();
    assemble_system_newmark(load, dt);
    x = llt.solve(f_x);
    y = llt.solve(f_y);
    z = llt.solve(f_z);
    update_multipliers();
    final_rel_update = compute_relative_state_update(
      x_iter_prev, x, y_iter_prev, y, z_iter_prev, z);
    final_max_pq_error = compute_max_pq_error();
    final_max_state_error = compute_max_state_update(
      x_iter_prev, x, y_iter_prev, y, z_iter_prev, z);

    if (final_rel_update < tol_outer) {
      converged = true;
      break;
    }
  }

  if (!converged) {
    ELFF_WARNING("EulerBeamDynamicInextensibleADDM::solve_newmark() final "
                 "relative state update = "
                 << final_rel_update << " at step " << time_iter << " after "
                 << iter << " iterations");
    ELFF_ABORT("EulerBeamDynamicInextensibleADDM::solve_newmark() did not "
               "converge.\n");
  }

  ELFF_LOG(time_iter << "\t" << final_max_pq_error << "\t"
                     << final_max_state_error << "\t" << iter);

  update_average_acceleration_state_component(x_old, x, vx_prev, ax_prev, dt);
  update_average_acceleration_state_component(y_old, y, vy_prev, ay_prev, dt);
  update_average_acceleration_state_component(z_old, z, vz_prev, az_prev, dt);
  apply_dynamic_state_boundary_conditions();

  x_prev = x;
  y_prev = y;
  z_prev = z;
  nodal_load_prev = load;
  have_prev_nodal_load = true;

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
  nodal_load_prev.clear();
  have_prev_uniform_load = false;
  have_prev_nodal_load = false;
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
  nodal_load_prev.clear();
  have_prev_uniform_load = false;
  have_prev_nodal_load = false;
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
  std::array<real_t, 3> load, real_t dt)
{
  clear_rhs();
  assemble_constraint_rhs();
  add_averaged_uniform_load_rhs(load);
  apply_midpoint_bending_rhs();

  const real_t coeff = 2.0 / (dt * dt);
  const VectorXd x_inertia = coeff * x_prev + (2.0 / dt) * vx_prev;
  const VectorXd y_inertia = coeff * y_prev + (2.0 / dt) * vy_prev;
  const VectorXd z_inertia = coeff * z_prev + (2.0 / dt) * vz_prev;

  f_x.noalias() += mass * x_inertia;
  f_y.noalias() += mass * y_inertia;
  f_z.noalias() += mass * z_inertia;

  A_unconstrained =
    coeff * mass + 0.5 * EI * K_bending + r_penalty * K_constraint;
  A = A_unconstrained;

  apply_boundary_condition_A();
  apply_boundary_condition_f();
  decompose_A();
}

void
EulerBeamDynamicInextensibleADDM::assemble_system_newmark(
  const std::vector<std::array<real_t, 3>>& load, real_t dt)
{
  clear_rhs();
  assemble_constraint_rhs();
  add_averaged_nodal_load_rhs(load);
  apply_midpoint_bending_rhs();

  const real_t coeff = 2.0 / (dt * dt);
  const VectorXd x_inertia = coeff * x_prev + (2.0 / dt) * vx_prev;
  const VectorXd y_inertia = coeff * y_prev + (2.0 / dt) * vy_prev;
  const VectorXd z_inertia = coeff * z_prev + (2.0 / dt) * vz_prev;

  f_x.noalias() += mass * x_inertia;
  f_y.noalias() += mass * y_inertia;
  f_z.noalias() += mass * z_inertia;

  A_unconstrained =
    coeff * mass + 0.5 * EI * K_bending + r_penalty * K_constraint;
  A = A_unconstrained;

  apply_boundary_condition_A();
  apply_boundary_condition_f();
  decompose_A();
}

void
EulerBeamDynamicInextensibleADDM::add_averaged_uniform_load_rhs(
  std::array<real_t, 3> load)
{
  std::array<real_t, 3> averaged = {
    0.5 * (load_prev[0] + load[0]),
    0.5 * (load_prev[1] + load[1]),
    0.5 * (load_prev[2] + load[2]),
  };
  add_uniform_load_rhs(averaged);
  add_point_boundary_loads();
}

void
EulerBeamDynamicInextensibleADDM::add_averaged_nodal_load_rhs(
  const std::vector<std::array<real_t, 3>>& load)
{
  ELFF_ASSERT(load.size() == mesh.get_nodes(),
              "Size of load vector must equal number of nodes.");
  ELFF_ASSERT(nodal_load_prev.size() == load.size(),
              "Stored nodal load history must match node count.");

  std::vector<std::array<real_t, 3>> averaged(load.size());
  for (size_t i = 0; i < load.size(); ++i) {
    averaged[i][0] = 0.5 * (nodal_load_prev[i][0] + load[i][0]);
    averaged[i][1] = 0.5 * (nodal_load_prev[i][1] + load[i][1]);
    averaged[i][2] = 0.5 * (nodal_load_prev[i][2] + load[i][2]);
  }
  add_nodal_load_rhs(averaged);
  add_point_boundary_loads();
}

void
EulerBeamDynamicInextensibleADDM::apply_midpoint_bending_rhs()
{
  f_x.noalias() -= 0.5 * EI * K_bending * x_prev;
  f_y.noalias() -= 0.5 * EI * K_bending * y_prev;
  f_z.noalias() -= 0.5 * EI * K_bending * z_prev;
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
EulerBeamDynamicInextensibleADDM::compute_relative_state_update(
  const VectorXd& x_old_iter,
  const VectorXd& x_new_iter,
  const VectorXd& y_old_iter,
  const VectorXd& y_new_iter,
  const VectorXd& z_old_iter,
  const VectorXd& z_new_iter) const
{
  const VectorXd dx = x_new_iter - x_old_iter;
  const VectorXd dy = y_new_iter - y_old_iter;
  const VectorXd dz = z_new_iter - z_old_iter;

  const real_t numerator =
    std::sqrt(dx.squaredNorm() + dy.squaredNorm() + dz.squaredNorm());
  const real_t denominator = std::max<real_t>(
    1e-14,
    std::sqrt(x_old_iter.squaredNorm() + y_old_iter.squaredNorm() +
              z_old_iter.squaredNorm()));

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
EulerBeamDynamicInextensibleADDM::compute_max_state_update(
  const VectorXd& x_old_iter,
  const VectorXd& x_new_iter,
  const VectorXd& y_old_iter,
  const VectorXd& y_new_iter,
  const VectorXd& z_old_iter,
  const VectorXd& z_new_iter) const
{
  const real_t max_dx = (x_new_iter - x_old_iter).cwiseAbs().maxCoeff();
  const real_t max_dy = (y_new_iter - y_old_iter).cwiseAbs().maxCoeff();
  const real_t max_dz = (z_new_iter - z_old_iter).cwiseAbs().maxCoeff();
  return std::max({ max_dx, max_dy, max_dz });
}

void
EulerBeamDynamicInextensibleADDM::update_average_acceleration_state_component(
  const VectorXd& u_old,
  const VectorXd& u_new,
  VectorXd& v_hist,
  VectorXd& a_hist,
  real_t dt)
{
  const VectorXd v_old = v_hist;
  v_hist = (2.0 / dt) * (u_new - u_old) - v_old;
  a_hist = (v_hist - v_old) / dt;
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
