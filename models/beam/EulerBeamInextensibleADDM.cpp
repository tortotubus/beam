#include "elff/models/beam/EulerBeamInextensibleADDM.hpp"
#include "elff/models/beam/EulerBeam.hpp"

namespace ELFF {
namespace Models {

EulerBeamInextensibleADDM::EulerBeamInextensibleADDM(real_t length, real_t EI,
                                                     size_t nodes,
                                                     EulerBeamBCs bcs,
                                                     real_t r_penalty)
    : EulerBeam(length, EI, nodes, bcs), elements(nodes - 1), dimension(3),
      dof((2 * nodes)), r_penalty(r_penalty), alpha(),
      A(SparseMatrix<real_t>(dof, dof)),
      A_unconstrained(SparseMatrix<real_t>(dof, dof)),
      K_bending(SparseMatrix<real_t>(dof, dof)),
      K_constraint(SparseMatrix<real_t>(dof, dof)),
      A_dense_cache(MatrixXd::Zero(dof, dof)),
      bending_element_matrix(Matrix<real_t, 4, 4>::Zero()),
      constraint_rhs_element_matrix(Matrix<real_t, 4, 3>::Zero()),
      midpoint_dH({0.0, 0.0, 0.0, 0.0}), x(VectorXd::Zero(dof)),
      y(VectorXd::Zero(dof)), z(VectorXd::Zero(dof)), f_x(VectorXd::Zero(dof)),
      f_y(VectorXd::Zero(dof)), f_z(VectorXd::Zero(dof)), llt(),
      lambda_x(VectorXd::Zero(nodes + elements)),
      lambda_y(VectorXd::Zero(nodes + elements)),
      lambda_z(VectorXd::Zero(nodes + elements)),
      p(VectorXd::Ones(nodes + elements)), q(VectorXd::Zero(nodes + elements)),
      r(VectorXd::Zero(nodes + elements)), xp(VectorXd::Zero(nodes + elements)),
      yp(VectorXd::Zero(nodes + elements)),
      zp(VectorXd::Zero(nodes + elements)), max_outer(100000), min_outer(1),
      omega_outer(1.0), tol_outer(1e-13), x_prev(VectorXd::Zero(dof)),
      y_prev(VectorXd::Zero(dof)), z_prev(VectorXd::Zero(dof)),
      vx_prev(VectorXd::Zero(dof)), vy_prev(VectorXd::Zero(dof)),
      vz_prev(VectorXd::Zero(dof)), ax_prev(VectorXd::Zero(dof)),
      ay_prev(VectorXd::Zero(dof)), az_prev(VectorXd::Zero(dof)),
      mass_matrix(SparseMatrix<real_t>(dof, dof)), load_prev({0., 0., 0.}),
      nodal_load_prev(), have_prev_uniform_load(false),
      have_prev_nodal_load(false) {
  apply_initial_condition();
  initialize_quadrature_cache();
  assemble_A();
  apply_boundary_condition_A();
  decompose_A();
  assemble_mass_matrix();
  x_prev = x;
  y_prev = y;
  z_prev = z;
}

EulerBeamInextensibleADDM::EulerBeamInextensibleADDM(real_t length, real_t EI,
                                                     real_t mu, size_t nodes,
                                                     EulerBeamBCs bcs,
                                                     real_t r_penalty)
    : EulerBeam(length, EI, mu, nodes, bcs), elements(nodes - 1), dimension(3),
      dof((2 * nodes)), r_penalty(r_penalty), alpha(),
      A(SparseMatrix<real_t>(dof, dof)),
      A_unconstrained(SparseMatrix<real_t>(dof, dof)),
      K_bending(SparseMatrix<real_t>(dof, dof)),
      K_constraint(SparseMatrix<real_t>(dof, dof)),
      A_dense_cache(MatrixXd::Zero(dof, dof)),
      bending_element_matrix(Matrix<real_t, 4, 4>::Zero()),
      constraint_rhs_element_matrix(Matrix<real_t, 4, 3>::Zero()),
      midpoint_dH({0.0, 0.0, 0.0, 0.0}), x(VectorXd::Zero(dof)),
      y(VectorXd::Zero(dof)), z(VectorXd::Zero(dof)), f_x(VectorXd::Zero(dof)),
      f_y(VectorXd::Zero(dof)), f_z(VectorXd::Zero(dof)), llt(),
      lambda_x(VectorXd::Zero(nodes + elements)),
      lambda_y(VectorXd::Zero(nodes + elements)),
      lambda_z(VectorXd::Zero(nodes + elements)),
      p(VectorXd::Ones(nodes + elements)), q(VectorXd::Zero(nodes + elements)),
      r(VectorXd::Zero(nodes + elements)), xp(VectorXd::Zero(nodes + elements)),
      yp(VectorXd::Zero(nodes + elements)),
      zp(VectorXd::Zero(nodes + elements)), max_outer(20000), min_outer(1),
      omega_outer(1.0), tol_outer(1e-9), x_prev(VectorXd::Zero(dof)),
      y_prev(VectorXd::Zero(dof)), z_prev(VectorXd::Zero(dof)),
      vx_prev(VectorXd::Zero(dof)),
      vy_prev(VectorXd::Zero(dof)), vz_prev(VectorXd::Zero(dof)),
      ax_prev(VectorXd::Zero(dof)), ay_prev(VectorXd::Zero(dof)),
      az_prev(VectorXd::Zero(dof)),
      mass_matrix(SparseMatrix<real_t>(dof, dof)),
      load_prev({0., 0., 0.}), nodal_load_prev(), have_prev_uniform_load(false),
      have_prev_nodal_load(false) {
  apply_initial_condition(mesh);
  initialize_quadrature_cache();
  assemble_A();
  apply_boundary_condition_A();
  decompose_A();
  assemble_mass_matrix();
  x_prev = x;
  y_prev = y;
  z_prev = z;
}

EulerBeamInextensibleADDM::EulerBeamInextensibleADDM(
    real_t length, real_t EI, real_t mu, size_t nodes,
    EulerBeamTimeDependentBCs bcs, real_t r_penalty)
    : EulerBeam(length, EI, mu, nodes, bcs), elements(nodes - 1), dimension(3),
      dof((2 * nodes)), r_penalty(r_penalty), alpha(),
      A(SparseMatrix<real_t>(dof, dof)),
      A_unconstrained(SparseMatrix<real_t>(dof, dof)),
      K_bending(SparseMatrix<real_t>(dof, dof)),
      K_constraint(SparseMatrix<real_t>(dof, dof)),
      A_dense_cache(MatrixXd::Zero(dof, dof)),
      bending_element_matrix(Matrix<real_t, 4, 4>::Zero()),
      constraint_rhs_element_matrix(Matrix<real_t, 4, 3>::Zero()),
      midpoint_dH({0.0, 0.0, 0.0, 0.0}), x(VectorXd::Zero(dof)),
      y(VectorXd::Zero(dof)), z(VectorXd::Zero(dof)), f_x(VectorXd::Zero(dof)),
      f_y(VectorXd::Zero(dof)), f_z(VectorXd::Zero(dof)), llt(),
      lambda_x(VectorXd::Zero(nodes + elements)),
      lambda_y(VectorXd::Zero(nodes + elements)),
      lambda_z(VectorXd::Zero(nodes + elements)),
      p(VectorXd::Ones(nodes + elements)), q(VectorXd::Zero(nodes + elements)),
      r(VectorXd::Zero(nodes + elements)), xp(VectorXd::Zero(nodes + elements)),
      yp(VectorXd::Zero(nodes + elements)),
      zp(VectorXd::Zero(nodes + elements)), max_outer(20000), min_outer(1),
      omega_outer(1.0), tol_outer(1e-9), x_prev(VectorXd::Zero(dof)),
      y_prev(VectorXd::Zero(dof)), z_prev(VectorXd::Zero(dof)),
      vx_prev(VectorXd::Zero(dof)),
      vy_prev(VectorXd::Zero(dof)), vz_prev(VectorXd::Zero(dof)),
      ax_prev(VectorXd::Zero(dof)), ay_prev(VectorXd::Zero(dof)),
      az_prev(VectorXd::Zero(dof)),
      mass_matrix(SparseMatrix<real_t>(dof, dof)),
      load_prev({0., 0., 0.}), nodal_load_prev(), have_prev_uniform_load(false),
      have_prev_nodal_load(false) {
  // Seed the live boundary-condition view from the prescribed history before
  // we assemble operators or initialize constrained state.
  apply_time_dependent_boundary_values(0);
  apply_initial_condition();
  apply_prescribed_boundary_values(x, y, z);
  initialize_quadrature_cache();
  apply_initial_condition_pq();
  assemble_A();
  apply_boundary_condition_A();
  decompose_A();
  assemble_mass_matrix();
  apply_time_dependent_boundary_kinematics(0);
  x_prev = x;
  y_prev = y;
  z_prev = z;
  apply_prescribed_boundary_values(x_prev, y_prev, z_prev);
  update_mesh();
}

EulerBeamInextensibleADDM::~EulerBeamInextensibleADDM() = default;

void EulerBeamInextensibleADDM::collect_boundary_dofs(
    std::vector<size_t> &idx, std::vector<real_t> &xvals,
    std::vector<real_t> &yvals, std::vector<real_t> &zvals) const {
  const size_t nodes = mesh.get_nodes();

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

    switch (bctype) {
    case free_bc:
    case point_force_bc:
    case point_torque_bc:
      break;
    case simple_bc:
      idx.push_back(2 * ni + 0);
      xvals.push_back(bcvals.position[0]);
      yvals.push_back(bcvals.position[1]);
      zvals.push_back(bcvals.position[2]);
      break;
    case clamped_bc:
      idx.push_back(2 * ni + 0);
      xvals.push_back(bcvals.position[0]);
      yvals.push_back(bcvals.position[1]);
      zvals.push_back(bcvals.position[2]);
      idx.push_back(2 * ni + 1);
      xvals.push_back(bcvals.slope[0]);
      yvals.push_back(bcvals.slope[1]);
      zvals.push_back(bcvals.slope[2]);
      break;
    }
  }
}

void EulerBeamInextensibleADDM::solve() { solve({0., 0., 0.}); }

void EulerBeamInextensibleADDM::solve(std::array<real_t, 3> load) {
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
    update_xy(load);
    x = (1.0 - omega_outer) * x_iter_prev + omega_outer * x;
    y = (1.0 - omega_outer) * y_iter_prev + omega_outer * y;
    z = (1.0 - omega_outer) * z_iter_prev + omega_outer * z;
    update_multipliers();
    final_rel_update = compute_relative_state_update(
        x_iter_prev, x, y_iter_prev, y, z_iter_prev, z);
    final_max_pq_error = compute_max_pq_error();
    final_max_state_error = compute_max_state_update(
        x_iter_prev, x, y_iter_prev, y, z_iter_prev, z);
    if (is_converged(x_iter_prev, x, y_iter_prev, y, z_iter_prev, z)) {
      converged = true;
      break;
    }
  }

  if (!converged) {
    ELFF_WARNING("EulerBeamInextensibleADDM::solve() final relative "
                 "state update = "
                 << final_rel_update << " after " << iter << " iterations");
  }

  ELFF_LOG(final_max_pq_error << " " << final_max_state_error << " " << iter);

  update_mesh();
}

void EulerBeamInextensibleADDM::solve(real_t dt) { solve(dt, {0., 0., 0.}); }

void EulerBeamInextensibleADDM::solve(real_t dt, std::array<real_t, 3> load) {
  const real_t alpha = 0.0;
  const real_t gamma = 0.5 - alpha;
  const real_t beta = 0.25 * (1.0 - alpha) * (1.0 - alpha);
  solve_newmark(dt, load, beta, gamma);
}

void EulerBeamInextensibleADDM::solve(real_t dt,
                                      std::vector<std::array<real_t, 3>> load) {
  ELFF_ASSERT(load.size() == mesh.get_nodes(),
              "Size of load vector must equal number of nodes.");

  const real_t alpha = 0.0;
  const real_t gamma = 0.5 - alpha;
  const real_t beta = 0.25 * (1.0 - alpha) * (1.0 - alpha);
  solve_newmark(dt, load, beta, gamma);
}

void EulerBeamInextensibleADDM::solve_averaged_load(
    real_t dt, const std::vector<std::array<real_t, 3>> &averaged_load) {
  ELFF_ASSERT(averaged_load.size() == mesh.get_nodes(),
              "Size of load vector must equal number of nodes.");

  const real_t alpha = 0.0;
  const real_t gamma = 0.5 - alpha;
  const real_t beta = 0.25 * (1.0 - alpha) * (1.0 - alpha);
  solve_newmark_averaged_load(dt, averaged_load, beta, gamma);
}

void EulerBeamInextensibleADDM::apply_initial_condition() {
  apply_initial_condition_xy();
  apply_initial_condition_pq();
  x_prev = x;
  y_prev = y;
  z_prev = z;
  vx_prev.setZero();
  vy_prev.setZero();
  vz_prev.setZero();
  ax_prev.setZero();
  ay_prev.setZero();
  az_prev.setZero();
  load_prev = {0., 0., 0.};
  nodal_load_prev.clear();
  have_prev_uniform_load = false;
  have_prev_nodal_load = false;
}

void EulerBeamInextensibleADDM::apply_initial_condition(EulerBeamMesh &bmesh) {
  apply_initial_condition_xy(bmesh);
  apply_initial_condition_pq();
  const auto velocities = bmesh.get_centerline_velocity();
  const auto accelerations = bmesh.get_centerline_acceleration();
  const auto slope_velocities = bmesh.get_slope_velocity();
  const auto slope_accelerations = bmesh.get_slope_acceleration();
  x_prev = x;
  y_prev = y;
  z_prev = z;
  vx_prev.setZero();
  vy_prev.setZero();
  vz_prev.setZero();
  ax_prev.setZero();
  ay_prev.setZero();
  az_prev.setZero();
  for (size_t i = 0; i < mesh.get_nodes(); ++i) {
    vx_prev(2 * i + 0) = velocities[i][0];
    vx_prev(2 * i + 1) = slope_velocities[i][0];
    vy_prev(2 * i + 0) = velocities[i][1];
    vy_prev(2 * i + 1) = slope_velocities[i][1];
    vz_prev(2 * i + 0) = velocities[i][2];
    vz_prev(2 * i + 1) = slope_velocities[i][2];

    ax_prev(2 * i + 0) = accelerations[i][0];
    ax_prev(2 * i + 1) = slope_accelerations[i][0];
    ay_prev(2 * i + 0) = accelerations[i][1];
    ay_prev(2 * i + 1) = slope_accelerations[i][1];
    az_prev(2 * i + 0) = accelerations[i][2];
    az_prev(2 * i + 1) = slope_accelerations[i][2];
  }
  load_prev = {0., 0., 0.};
  nodal_load_prev.clear();
  have_prev_uniform_load = false;
  have_prev_nodal_load = false;
}

void EulerBeamInextensibleADDM::set_time_dependent_boundary_conditions(
    EulerBeamTimeDependentBCs bcs) {
  EulerBeam::set_time_dependent_boundary_conditions(bcs);
  apply_time_dependent_boundary_values(0);
  apply_prescribed_boundary_values(x, y, z);
  apply_prescribed_boundary_values(x_prev, y_prev, z_prev);
  apply_time_dependent_boundary_kinematics(0);
  apply_boundary_condition_lambda();
  apply_boundary_condition_A();
  decompose_A();
  update_mesh();
}

const VectorXd &EulerBeamInextensibleADDM::get_lambda_x() const {
  return lambda_x;
}

const VectorXd &EulerBeamInextensibleADDM::get_lambda_y() const {
  return lambda_y;
}

const VectorXd &EulerBeamInextensibleADDM::get_lambda_z() const {
  return lambda_z;
}

const MatrixXd &EulerBeamInextensibleADDM::get_A() const {
  A_dense_cache = MatrixXd(A);
  return A_dense_cache;
}

void EulerBeamInextensibleADDM::update_mesh() {
  const size_t nodes = this->mesh.get_nodes();
  std::vector<std::array<real_t, 3>> &centerline = mesh.get_centerline();
  std::vector<std::array<real_t, 3>> &slope = mesh.get_slope();
  std::vector<std::array<real_t, 3>> &velocity = mesh.get_centerline_velocity();
  std::vector<std::array<real_t, 3>> &acceleration =
      mesh.get_centerline_acceleration();
  std::vector<std::array<real_t, 3>> &slope_velocity =
      mesh.get_slope_velocity();
  std::vector<std::array<real_t, 3>> &slope_acceleration =
      mesh.get_slope_acceleration();
  std::vector<real_t> &s = mesh.get_curvilinear_axis();
  (void)s;

  for (size_t i = 0; i < nodes; ++i) {
    centerline[i][0] = x(2 * i);
    centerline[i][1] = y(2 * i);
    centerline[i][2] = z(2 * i);
    slope[i][0] = x(2 * i + 1);
    slope[i][1] = y(2 * i + 1);
    slope[i][2] = z(2 * i + 1);
    velocity[i][0] = vx_prev(2 * i + 0);
    velocity[i][1] = vy_prev(2 * i + 0);
    velocity[i][2] = vz_prev(2 * i + 0);
    acceleration[i][0] = ax_prev(2 * i + 0);
    acceleration[i][1] = ay_prev(2 * i + 0);
    acceleration[i][2] = az_prev(2 * i + 0);
    slope_velocity[i][0] = vx_prev(2 * i + 1);
    slope_velocity[i][1] = vy_prev(2 * i + 1);
    slope_velocity[i][2] = vz_prev(2 * i + 1);
    slope_acceleration[i][0] = ax_prev(2 * i + 1);
    slope_acceleration[i][1] = ay_prev(2 * i + 1);
    slope_acceleration[i][2] = az_prev(2 * i + 1);
  }
}

void EulerBeamInextensibleADDM::apply_time_dependent_boundary_values(
    size_t step_idx) {
  ELFF_ASSERT(have_time_dependent_boundary_conditions,
              "Time-dependent boundary conditions were not provided.");

  const size_t hist_idx =
      time_dependent_boundary_conditions.time_zero_idx + step_idx;

  for (size_t bi = 0; bi < 2; ++bi) {
    boundary_conditions.end[bi] = time_dependent_boundary_conditions.end[bi];
    boundary_conditions.type[bi] = time_dependent_boundary_conditions.type[bi];

    const EulerBeamBCType bctype = boundary_conditions.type[bi];
    const auto &history = time_dependent_boundary_conditions.history[bi];
    EulerBeamBCVals &vals = boundary_conditions.vals[bi];

    switch (bctype) {
    case free_bc:
      break;
    case simple_bc:
      ELFF_ASSERT(hist_idx < history.position_history.size(),
                  "Not enough prescribed boundary position history for ADDM.");
      for (size_t d = 0; d < 3; ++d) {
        vals.position[d] = history.position_history[hist_idx][d];
      }
      break;
    case clamped_bc:
      ELFF_ASSERT(hist_idx < history.position_history.size(),
                  "Not enough prescribed boundary position history for ADDM.");
      ELFF_ASSERT(hist_idx < history.slope_history.size(),
                  "Not enough prescribed boundary slope history for ADDM.");
      for (size_t d = 0; d < 3; ++d) {
        vals.position[d] = history.position_history[hist_idx][d];
        vals.slope[d] = history.slope_history[hist_idx][d];
      }
      break;
    case point_force_bc:
      ELFF_ASSERT(hist_idx < history.force_history.size(),
                  "Not enough prescribed boundary force history for ADDM.");
      for (size_t d = 0; d < 3; ++d) {
        vals.force[d] = history.force_history[hist_idx][d];
      }
      break;
    case point_torque_bc:
      ELFF_ASSERT(hist_idx < history.torque_history.size(),
                  "Not enough prescribed boundary torque history for ADDM.");
      for (size_t d = 0; d < 3; ++d) {
        vals.torque[d] = history.torque_history[hist_idx][d];
      }
      break;
    }
  }
}

void EulerBeamInextensibleADDM::apply_time_dependent_boundary_kinematics(
    size_t step_idx) {
  if (!have_time_dependent_boundary_conditions) {
    clear_constrained_dynamic_history();
    return;
  }

  const size_t hist_idx =
      time_dependent_boundary_conditions.time_zero_idx + step_idx;
  const size_t nodes = mesh.get_nodes();

  for (size_t bi = 0; bi < 2; ++bi) {
    const EulerBeamBCType bctype = boundary_conditions.type[bi];
    if (bctype != simple_bc && bctype != clamped_bc) {
      continue;
    }

    const auto &history = time_dependent_boundary_conditions.history[bi];
    ELFF_ASSERT(hist_idx < history.velocity_history.size(),
                "Not enough prescribed boundary velocity history for ADDM.");
    ELFF_ASSERT(
        hist_idx < history.acceleration_history.size(),
        "Not enough prescribed boundary acceleration history for ADDM.");

    const size_t ni = (boundary_conditions.end[bi] == left) ? 0 : nodes - 1;
    vx_prev(2 * ni + 0) = history.velocity_history[hist_idx][0];
    vy_prev(2 * ni + 0) = history.velocity_history[hist_idx][1];
    vz_prev(2 * ni + 0) = history.velocity_history[hist_idx][2];
    ax_prev(2 * ni + 0) = history.acceleration_history[hist_idx][0];
    ay_prev(2 * ni + 0) = history.acceleration_history[hist_idx][1];
    az_prev(2 * ni + 0) = history.acceleration_history[hist_idx][2];

    if (bctype == clamped_bc) {
      ELFF_ASSERT(hist_idx < history.slope_velocity_history.size(),
                  "Not enough prescribed boundary slope velocity history for "
                  "ADDM.");
      ELFF_ASSERT(
          hist_idx < history.slope_acceleration_history.size(),
          "Not enough prescribed boundary slope acceleration history for "
          "ADDM.");

      vx_prev(2 * ni + 1) = history.slope_velocity_history[hist_idx][0];
      vy_prev(2 * ni + 1) = history.slope_velocity_history[hist_idx][1];
      vz_prev(2 * ni + 1) = history.slope_velocity_history[hist_idx][2];
      ax_prev(2 * ni + 1) = history.slope_acceleration_history[hist_idx][0];
      ay_prev(2 * ni + 1) = history.slope_acceleration_history[hist_idx][1];
      az_prev(2 * ni + 1) = history.slope_acceleration_history[hist_idx][2];
    }
  }
}

void EulerBeamInextensibleADDM::apply_prescribed_boundary_values(
    VectorXd &x_state, VectorXd &y_state, VectorXd &z_state) const {
  std::vector<size_t> idx;
  std::vector<real_t> xvals, yvals, zvals;
  collect_boundary_dofs(idx, xvals, yvals, zvals);

  for (size_t i = 0; i < idx.size(); ++i) {
    const size_t d = idx[i];
    x_state(d) = xvals[i];
    y_state(d) = yvals[i];
    z_state(d) = zvals[i];
  }
}

void EulerBeamInextensibleADDM::clear_constrained_dynamic_history() {
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

void EulerBeamInextensibleADDM::apply_initial_condition_xy() {
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

void EulerBeamInextensibleADDM::apply_initial_condition_xy(
    EulerBeamMesh &bmesh) {
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

void EulerBeamInextensibleADDM::initialize_quadrature_cache() {
  const real_t h = mesh.get_ds();
  const real_t xi_q[3] = {0.1127016654, 0.5, 0.8872983346};
  const real_t w_q[3] = {0.2777777778, 0.4444444444, 0.2777777778};

  midpoint_dH = ELFF::FEM::CubicHermite<real_t>::derivs(0.5, h);
  constraint_rhs_element_matrix.setZero();

  for (size_t qi = 0; qi < 3; ++qi) {
    const auto L = ELFF::FEM::QuadraticLagrange<real_t>::values(xi_q[qi]);
    const auto dH = ELFF::FEM::CubicHermite<real_t>::derivs(xi_q[qi], h);
    const real_t weight = w_q[qi] * h;

    for (size_t a = 0; a < 4; ++a) {
      for (size_t j = 0; j < 3; ++j) {
        constraint_rhs_element_matrix(a, j) += L[j] * dH[a] * weight;
      }
    }
  }
}

void EulerBeamInextensibleADDM::compute_slopes_collocation() {
  xp.setZero();
  yp.setZero();
  zp.setZero();
  const size_t nodes = this->mesh.get_nodes();

  for (size_t ni = 0; ni < nodes; ni++) {
    xp[ni] = x[2 * ni + 1];
    yp[ni] = y[2 * ni + 1];
    zp[ni] = z[2 * ni + 1];
  }

  for (size_t ei = 0; ei < elements; ei++) {
    const size_t edofs[4] = {2 * (ei + 0) + 0, 2 * (ei + 0) + 1,
                             2 * (ei + 1) + 0, 2 * (ei + 1) + 1};

    xp[nodes + ei] =
        midpoint_dH[0] * x[edofs[0]] + midpoint_dH[1] * x[edofs[1]] +
        midpoint_dH[2] * x[edofs[2]] + midpoint_dH[3] * x[edofs[3]];
    yp[nodes + ei] =
        midpoint_dH[0] * y[edofs[0]] + midpoint_dH[1] * y[edofs[1]] +
        midpoint_dH[2] * y[edofs[2]] + midpoint_dH[3] * y[edofs[3]];
    zp[nodes + ei] =
        midpoint_dH[0] * z[edofs[0]] + midpoint_dH[1] * z[edofs[1]] +
        midpoint_dH[2] * z[edofs[2]] + midpoint_dH[3] * z[edofs[3]];
  }
}

void EulerBeamInextensibleADDM::apply_initial_condition_pq() {
  const size_t nodes = this->mesh.get_nodes();
  compute_slopes_collocation();
  for (size_t ci = 0; ci < nodes + elements; ci++) {
    switch (dimension) {
    case 2: {
      const real_t xprime = xp[ci];
      const real_t yprime = yp[ci];
      const real_t den =
          std::max(1e-14, sqrt(xprime * xprime + yprime * yprime));
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
}

void EulerBeamInextensibleADDM::apply_boundary_condition_lambda() {
  const size_t nodes = mesh.get_nodes();
  for (size_t bi = 0; bi < 2; ++bi) {
    const EulerBeamBCType bctype = boundary_conditions.type[bi];
    if (bctype == free_bc) {
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
}

void EulerBeamInextensibleADDM::update_pq() {
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
}

void EulerBeamInextensibleADDM::apply_boundary_condition_A() {
  A = A_unconstrained;

  std::vector<size_t> idx;
  std::vector<real_t> xvals, yvals, zvals;
  collect_boundary_dofs(idx, xvals, yvals, zvals);

  std::vector<char> constrained(dof, 0);
  for (size_t d : idx) {
    constrained[d] = 1;
  }

  // Zero all couplings touching constrained rows/cols in one sparse sweep.
  for (int col = 0; col < A.outerSize(); ++col) {
    for (SparseMatrix<real_t>::InnerIterator it(A, col); it; ++it) {
      if ((constrained[it.row()] || constrained[it.col()]) &&
          it.row() != it.col()) {
        it.valueRef() = 0.0;
      }
    }
  }

  for (size_t d : idx) {
    A.coeffRef(d, d) = 1.0;
  }

  A.prune(0.0);
  A.makeCompressed();
}

void EulerBeamInextensibleADDM::assemble_A() {
  using Tpl = Triplet<real_t>;

  const real_t h = this->mesh.get_ds();
  const real_t xi_q[3] = {0.1127016654, 0.5, 0.8872983346};
  const real_t w_q[3] = {0.2777777778, 0.4444444444, 0.2777777778};

  bending_element_matrix.setZero();
  Matrix<real_t, 4, 4> constraint_element_matrix = Matrix<real_t, 4, 4>::Zero();

  for (size_t q = 0; q < 3; ++q) {
    const real_t xi = xi_q[q];
    const real_t w = w_q[q];

    const auto dH = ELFF::FEM::CubicHermite<real_t>::derivs(xi, h);
    const auto ddH = ELFF::FEM::CubicHermite<real_t>::second_derivs(xi, h);

    for (size_t a = 0; a < 4; ++a) {
      for (size_t b = 0; b < 4; ++b) {
        bending_element_matrix(a, b) += ddH[a] * ddH[b] * w * h;
        constraint_element_matrix(a, b) += dH[a] * dH[b] * w * h;
      }
    }
  }

  std::vector<Tpl> bending_triplets;
  std::vector<Tpl> constraint_triplets;
  bending_triplets.reserve(elements * 16);
  constraint_triplets.reserve(elements * 16);

  for (size_t e = 0; e < this->elements; ++e) {
    const size_t edofs[4] = {2 * (e + 0) + 0, 2 * (e + 0) + 1, 2 * (e + 1) + 0,
                             2 * (e + 1) + 1};

    for (size_t a = 0; a < 4; ++a) {
      for (size_t b = 0; b < 4; ++b) {
        bending_triplets.emplace_back(edofs[a], edofs[b],
                                      bending_element_matrix(a, b));
        constraint_triplets.emplace_back(edofs[a], edofs[b],
                                         constraint_element_matrix(a, b));
      }
    }
  }

  K_bending.resize(dof, dof);
  K_constraint.resize(dof, dof);
  K_bending.setFromTriplets(bending_triplets.begin(), bending_triplets.end(),
                            std::plus<real_t>());
  K_constraint.setFromTriplets(constraint_triplets.begin(),
                               constraint_triplets.end(), std::plus<real_t>());
  K_bending.makeCompressed();
  K_constraint.makeCompressed();

  this->A_unconstrained = this->EI * K_bending + this->r_penalty * K_constraint;
  this->A_unconstrained.makeCompressed();
  this->A = this->A_unconstrained;
}

void EulerBeamInextensibleADDM::decompose_A() {
  llt.compute(A);
  ELFF_ASSERT(llt.info() == Success,
              "Sparse factorization failed for EulerBeamInextensibleADDM.\n");
}

void EulerBeamInextensibleADDM::clear_rhs() {
  f_x.setZero();
  f_y.setZero();
  f_z.setZero();
}

void EulerBeamInextensibleADDM::assemble_constraint_rhs() {
  const size_t nodes = this->mesh.get_nodes();

  for (size_t e = 0; e < elements; ++e) {
    const size_t edofs[4] = {2 * (e + 0) + 0, 2 * (e + 0) + 1, 2 * (e + 1) + 0,
                             2 * (e + 1) + 1};
    const size_t li = e;
    const size_t mi = nodes + e;
    const size_t ri = e + 1;

    Matrix<real_t, 3, 1> cx_elem;
    Matrix<real_t, 3, 1> cy_elem;
    Matrix<real_t, 3, 1> cz_elem;
    cx_elem << lambda_x[li] + r_penalty * p[li],
        lambda_x[mi] + r_penalty * p[mi], lambda_x[ri] + r_penalty * p[ri];
    cy_elem << lambda_y[li] + r_penalty * q[li],
        lambda_y[mi] + r_penalty * q[mi], lambda_y[ri] + r_penalty * q[ri];
    cz_elem << lambda_z[li] + r_penalty * r[li],
        lambda_z[mi] + r_penalty * r[mi], lambda_z[ri] + r_penalty * r[ri];

    const Matrix<real_t, 4, 1> fxe = constraint_rhs_element_matrix * cx_elem;
    const Matrix<real_t, 4, 1> fye = constraint_rhs_element_matrix * cy_elem;
    const Matrix<real_t, 4, 1> fze = constraint_rhs_element_matrix * cz_elem;

    for (size_t a = 0; a < 4; ++a) {
      f_x(edofs[a]) += fxe(a);
      f_y(edofs[a]) += fye(a);
      f_z(edofs[a]) += fze(a);
    }
  }
}

void EulerBeamInextensibleADDM::add_uniform_load_rhs(
    std::array<real_t, 3> load) {
  const real_t h = this->mesh.get_ds();

  const real_t xi_q[3] = {0.1127016654, 0.5, 0.8872983346};
  const real_t w_q[3] = {0.2777777778, 0.4444444444, 0.2777777778};

  for (size_t e = 0; e < elements; ++e) {
    const size_t edofs[4] = {2 * (e + 0) + 0, 2 * (e + 0) + 1, 2 * (e + 1) + 0,
                             2 * (e + 1) + 1};
    real_t fxe[4] = {0, 0, 0, 0};
    real_t fye[4] = {0, 0, 0, 0};
    real_t fze[4] = {0, 0, 0, 0};

    for (size_t qi = 0; qi < 3; ++qi) {
      const real_t xi = xi_q[qi];
      const real_t w = w_q[qi];

      const auto H = ELFF::FEM::CubicHermite<real_t>::values(xi, h);

      for (size_t a = 0; a < 4; ++a) {
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

void EulerBeamInextensibleADDM::add_nodal_load_rhs(
    const std::vector<std::array<real_t, 3>> &load) {
  ELFF_ASSERT(load.size() == mesh.get_nodes(),
              "Size of load vector must equal number of nodes.");

  const real_t h = this->mesh.get_ds();

  const real_t xi_q[3] = {0.1127016654, 0.5, 0.8872983346};
  const real_t w_q[3] = {0.2777777778, 0.4444444444, 0.2777777778};

  for (size_t e = 0; e < elements; ++e) {
    const size_t edofs[4] = {2 * (e + 0) + 0, 2 * (e + 0) + 1, 2 * (e + 1) + 0,
                             2 * (e + 1) + 1};
    real_t fxe[4] = {0, 0, 0, 0};
    real_t fye[4] = {0, 0, 0, 0};
    real_t fze[4] = {0, 0, 0, 0};

    for (size_t qi = 0; qi < 3; ++qi) {
      const real_t xi = xi_q[qi];
      const real_t w = w_q[qi];

      const auto M = ELFF::FEM::LinearShape<real_t>::values(xi);
      const auto H = ELFF::FEM::CubicHermite<real_t>::values(xi, h);

      const real_t fx_q = M[0] * load[e][0] + M[1] * load[e + 1][0];
      const real_t fy_q = M[0] * load[e][1] + M[1] * load[e + 1][1];
      const real_t fz_q = M[0] * load[e][2] + M[1] * load[e + 1][2];

      for (size_t a = 0; a < 4; ++a) {
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
}

void EulerBeamInextensibleADDM::add_point_boundary_loads() {
  const size_t nodes_count = mesh.get_nodes();

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
      ni = nodes_count - 1;
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

void EulerBeamInextensibleADDM::assemble_f(std::array<real_t, 3> load) {
  clear_rhs();
  assemble_constraint_rhs();
  add_uniform_load_rhs(load);
  add_point_boundary_loads();
}

void EulerBeamInextensibleADDM::apply_boundary_condition_f() {
  std::vector<size_t> idx;
  std::vector<real_t> xvals, yvals, zvals;
  collect_boundary_dofs(idx, xvals, yvals, zvals);

  // Reproduce the original sequential elimination exactly, but read each
  // constrained column from the sparse unconstrained operator and ignore rows
  // that were already eliminated by earlier constrained DOFs.
  std::vector<char> eliminated(dof, 0);

  for (size_t i = 0; i < idx.size(); i++) {
    const int d = static_cast<int>(idx[i]);
    for (SparseMatrix<real_t>::InnerIterator it(A_unconstrained, d); it; ++it) {
      const int row = it.row();
      if (eliminated[row]) {
        continue;
      }

      const real_t value = it.value();
      f_x(row) -= value * xvals[i];
      f_y(row) -= value * yvals[i];
      f_z(row) -= value * zvals[i];
    }

    f_x(d) = xvals[i];
    f_y(d) = yvals[i];
    f_z(d) = zvals[i];
    eliminated[d] = 1;
  }
}

void EulerBeamInextensibleADDM::update_xy(std::array<real_t, 3> load) {
  assemble_f(load);
  apply_boundary_condition_f();

  x = llt.solve(f_x);
  y = llt.solve(f_y);
  z = llt.solve(f_z);
}

void EulerBeamInextensibleADDM::update_multipliers() {
  const size_t nodes = mesh.get_nodes();
  compute_slopes_collocation();

  for (size_t ci = 0; ci < nodes + elements; ci++) {
    lambda_x[ci] += r_penalty * (p[ci] - xp[ci]);
    lambda_y[ci] += r_penalty * (q[ci] - yp[ci]);
    lambda_z[ci] += r_penalty * (r[ci] - zp[ci]);
  }

  apply_boundary_condition_lambda();
}

bool EulerBeamInextensibleADDM::is_converged(const VectorXd &x_old_iter,
                                             const VectorXd &x_new_iter,
                                             const VectorXd &y_old_iter,
                                             const VectorXd &y_new_iter,
                                             const VectorXd &z_old_iter,
                                             const VectorXd &z_new_iter) const {
  return compute_relative_state_update(x_old_iter, x_new_iter, y_old_iter,
                                       y_new_iter, z_old_iter,
                                       z_new_iter) < tol_outer;
}

void EulerBeamInextensibleADDM::solve_newmark(real_t dt,
                                              std::array<real_t, 3> load,
                                              real_t beta, real_t gamma) {
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
              "EulerBeamInextensibleADDM implements the average-"
              "acceleration Newmark variant used in the paper.\n");

  if (!have_prev_uniform_load) {
    load_prev = load;
    have_prev_uniform_load = true;
  }
  have_prev_nodal_load = false;

  x = x_prev;
  y = y_prev;
  z = z_prev;
  if (have_time_dependent_boundary_conditions) {
    apply_time_dependent_boundary_values(time_iter + 1);
    apply_prescribed_boundary_values(x, y, z);
  }

  const VectorXd x_old = x_prev;
  const VectorXd y_old = y_prev;
  const VectorXd z_old = z_prev;
  prepare_system_newmark(dt);
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
    assemble_system_newmark_rhs(load, dt);
    const VectorXd x_raw = llt.solve(f_x);
    const VectorXd y_raw = llt.solve(f_y);
    const VectorXd z_raw = llt.solve(f_z);
    x = (1.0 - omega_outer) * x_iter_prev + omega_outer * x_raw;
    y = (1.0 - omega_outer) * y_iter_prev + omega_outer * y_raw;
    z = (1.0 - omega_outer) * z_iter_prev + omega_outer * z_raw;
    update_multipliers();
    final_rel_update = compute_relative_state_update(
        x_iter_prev, x, y_iter_prev, y, z_iter_prev, z);
    final_max_pq_error = compute_max_pq_error();
    final_max_state_error = compute_max_state_update(
        x_iter_prev, x, y_iter_prev, y, z_iter_prev, z);

    if (iter >= min_outer) {
      if (is_converged(x_iter_prev, x, y_iter_prev, y, z_iter_prev, z)) {
        converged = true;
        break;
      }
    }
  }

  if (!converged) {
    ELFF_WARNING("EulerBeamInextensibleADDM::solve_newmark() final "
                 "relative state update = "
                 << final_rel_update << " at step " << time_iter << " after "
                 << iter << " iterations");
  }

  ELFF_LOG(time_iter << " " << final_rel_update << " " << final_max_state_error
                     << " " << iter);

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

void EulerBeamInextensibleADDM::solve_newmark(
    real_t dt, std::vector<std::array<real_t, 3>> load, real_t beta,
    real_t gamma) {
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
              "EulerBeamInextensibleADDM implements the average-"
              "acceleration Newmark variant used in the paper.\n");

  if (!have_prev_nodal_load) {
    nodal_load_prev = load;
    have_prev_nodal_load = true;
  }
  have_prev_uniform_load = false;

  x = x_prev;
  y = y_prev;
  z = z_prev;
  if (have_time_dependent_boundary_conditions) {
    apply_time_dependent_boundary_values(time_iter + 1);
    apply_prescribed_boundary_values(x, y, z);
  }

  const VectorXd x_old = x_prev;
  const VectorXd y_old = y_prev;
  const VectorXd z_old = z_prev;
  prepare_system_newmark(dt);
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
    assemble_system_newmark_rhs(load, dt);
    const VectorXd x_raw = llt.solve(f_x);
    const VectorXd y_raw = llt.solve(f_y);
    const VectorXd z_raw = llt.solve(f_z);
    x = (1.0 - omega_outer) * x_iter_prev + omega_outer * x_raw;
    y = (1.0 - omega_outer) * y_iter_prev + omega_outer * y_raw;
    z = (1.0 - omega_outer) * z_iter_prev + omega_outer * z_raw;
    update_multipliers();
    final_rel_update = compute_relative_state_update(
        x_iter_prev, x, y_iter_prev, y, z_iter_prev, z);
    final_max_pq_error = compute_max_pq_error();
    final_max_state_error = compute_max_state_update(
        x_iter_prev, x, y_iter_prev, y, z_iter_prev, z);

    if (iter >= min_outer) {
      if (is_converged(x_iter_prev, x, y_iter_prev, y, z_iter_prev, z)) {
        converged = true;
        break;
      }
    }
  }

  if (!converged) {
    ELFF_WARNING("EulerBeamInextensibleADDM::solve_newmark() final "
                 "relative state update = "
                 << final_rel_update << " at step " << time_iter << " after "
                 << iter << " iterations");
    // ELFF_ABORT("EulerBeamInextensibleADDM::solve_newmark() did not "
    //            "converge.\n");
  }

  ELFF_LOG(time_iter << " " << final_rel_update << " " << final_max_state_error
                     << " " << iter);

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

void EulerBeamInextensibleADDM::solve_newmark_averaged_load(
    real_t dt, const std::vector<std::array<real_t, 3>> &averaged_load,
    real_t beta, real_t gamma) {
  ELFF_ASSERT(averaged_load.size() == mesh.get_nodes(),
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
              "EulerBeamInextensibleADDM implements the average-"
              "acceleration Newmark variant used in the paper.\n");

  have_prev_uniform_load = false;
  have_prev_nodal_load = false;

  x = x_prev;
  y = y_prev;
  z = z_prev;
  if (have_time_dependent_boundary_conditions) {
    apply_time_dependent_boundary_values(time_iter + 1);
    apply_prescribed_boundary_values(x, y, z);
  }

  const VectorXd x_old = x_prev;
  const VectorXd y_old = y_prev;
  const VectorXd z_old = z_prev;
  prepare_system_newmark(dt);
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
    assemble_system_newmark_rhs_averaged(averaged_load, dt);
    const VectorXd x_raw = llt.solve(f_x);
    const VectorXd y_raw = llt.solve(f_y);
    const VectorXd z_raw = llt.solve(f_z);
    x = (1.0 - omega_outer) * x_iter_prev + omega_outer * x_raw;
    y = (1.0 - omega_outer) * y_iter_prev + omega_outer * y_raw;
    z = (1.0 - omega_outer) * z_iter_prev + omega_outer * z_raw;
    update_multipliers();
    final_rel_update = compute_relative_state_update(
        x_iter_prev, x, y_iter_prev, y, z_iter_prev, z);
    final_max_pq_error = compute_max_pq_error();
    final_max_state_error = compute_max_state_update(
        x_iter_prev, x, y_iter_prev, y, z_iter_prev, z);

    if (iter >= min_outer) {
      if (is_converged(x_iter_prev, x, y_iter_prev, y, z_iter_prev, z)) {
        converged = true;
        break;
      }
    }
  }

  if (!converged) {
    ELFF_WARNING("EulerBeamInextensibleADDM::solve_newmark_averaged_load() "
                 "final relative state update = "
                 << final_rel_update << " at step " << time_iter << " after "
                 << iter << " iterations");
  }

  ELFF_LOG(time_iter << " " << final_rel_update << " " << final_max_state_error
                     << " " << iter);

  update_average_acceleration_state_component(x_old, x, vx_prev, ax_prev, dt);
  update_average_acceleration_state_component(y_old, y, vy_prev, ay_prev, dt);
  update_average_acceleration_state_component(z_old, z, vz_prev, az_prev, dt);
  apply_dynamic_state_boundary_conditions();

  x_prev = x;
  y_prev = y;
  z_prev = z;
  nodal_load_prev.clear();

  update_mesh();

  ++time_iter;
  t += dt;
}

void EulerBeamInextensibleADDM::assemble_mass_matrix() {
  using Tpl = Triplet<real_t>;

  std::vector<Tpl> triplets;
  triplets.reserve(elements * 16);

  const real_t h = mesh.get_ds();
  const real_t xi_q[3] = {0.1127016654, 0.5, 0.8872983346};
  const real_t w_q[3] = {0.2777777778, 0.4444444444, 0.2777777778};

  for (size_t e = 0; e < elements; ++e) {
    const std::array<size_t, 4> edofs = {2 * (e + 0) + 0, 2 * (e + 0) + 1,
                                         2 * (e + 1) + 0, 2 * (e + 1) + 1};

    real_t me[4][4] = {{0.0}};

    for (size_t qi = 0; qi < 3; ++qi) {
      const real_t xi = xi_q[qi];
      const real_t w = w_q[qi];
      const auto H = ELFF::FEM::CubicHermite<real_t>::values(xi, h);

      for (size_t a = 0; a < 4; ++a) {
        for (size_t b = 0; b < 4; ++b) {
          me[a][b] += mu * H[a] * H[b] * w * h;
        }
      }
    }

    for (size_t a = 0; a < 4; ++a) {
      for (size_t b = 0; b < 4; ++b) {
        if (me[a][b] == 0.0) {
          continue;
        }
        triplets.emplace_back(edofs[a], edofs[b], me[a][b]);
      }
    }
  }

  mass_matrix.resize(dof, dof);
  mass_matrix.setFromTriplets(triplets.begin(), triplets.end());
  mass_matrix.makeCompressed();
}

void EulerBeamInextensibleADDM::prepare_system_newmark(real_t dt) {
  const real_t coeff = 2.0 / (dt * dt);

  A_unconstrained = 0.5 * EI * K_bending + r_penalty * K_constraint
                  + coeff * mass_matrix;
  A_unconstrained.makeCompressed();
  A = A_unconstrained;

  apply_boundary_condition_A();
  decompose_A();
}

void EulerBeamInextensibleADDM::assemble_system_newmark_rhs(
    std::array<real_t, 3> load, real_t dt) {
  clear_rhs();
  assemble_constraint_rhs();
  add_averaged_uniform_load_rhs(load);
  apply_midpoint_bending_rhs();

  const real_t coeff = 2.0 / (dt * dt);
  const VectorXd x_inertia = coeff * x_prev + (2.0 / dt) * vx_prev;
  const VectorXd y_inertia = coeff * y_prev + (2.0 / dt) * vy_prev;
  const VectorXd z_inertia = coeff * z_prev + (2.0 / dt) * vz_prev;

  f_x += mass_matrix * x_inertia;
  f_y += mass_matrix * y_inertia;
  f_z += mass_matrix * z_inertia;

  apply_boundary_condition_f();
}

void EulerBeamInextensibleADDM::assemble_system_newmark_rhs(
    const std::vector<std::array<real_t, 3>> &load, real_t dt) {
  clear_rhs();
  assemble_constraint_rhs();
  add_averaged_nodal_load_rhs(load);
  apply_midpoint_bending_rhs();

  const real_t coeff = 2.0 / (dt * dt);
  const VectorXd x_inertia = coeff * x_prev + (2.0 / dt) * vx_prev;
  const VectorXd y_inertia = coeff * y_prev + (2.0 / dt) * vy_prev;
  const VectorXd z_inertia = coeff * z_prev + (2.0 / dt) * vz_prev;

  f_x += mass_matrix * x_inertia;
  f_y += mass_matrix * y_inertia;
  f_z += mass_matrix * z_inertia;

  apply_boundary_condition_f();
}

void EulerBeamInextensibleADDM::assemble_system_newmark_rhs_averaged(
    const std::vector<std::array<real_t, 3>> &averaged_load, real_t dt) {
  clear_rhs();
  assemble_constraint_rhs();
  add_nodal_load_rhs(averaged_load);
  add_point_boundary_loads();
  apply_midpoint_bending_rhs();

  const real_t coeff = 2.0 / (dt * dt);
  const VectorXd x_inertia = coeff * x_prev + (2.0 / dt) * vx_prev;
  const VectorXd y_inertia = coeff * y_prev + (2.0 / dt) * vy_prev;
  const VectorXd z_inertia = coeff * z_prev + (2.0 / dt) * vz_prev;

  f_x += mass_matrix * x_inertia;
  f_y += mass_matrix * y_inertia;
  f_z += mass_matrix * z_inertia;

  apply_boundary_condition_f();
}

void EulerBeamInextensibleADDM::add_averaged_uniform_load_rhs(
    std::array<real_t, 3> load) {
  std::array<real_t, 3> averaged = {
      0.5 * (load_prev[0] + load[0]),
      0.5 * (load_prev[1] + load[1]),
      0.5 * (load_prev[2] + load[2]),
  };
  add_uniform_load_rhs(averaged);
  add_point_boundary_loads();
}

void EulerBeamInextensibleADDM::add_averaged_nodal_load_rhs(
    const std::vector<std::array<real_t, 3>> &load) {
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

void EulerBeamInextensibleADDM::apply_midpoint_bending_rhs() {
  const real_t coeff = -0.5 * EI;

  for (size_t e = 0; e < elements; ++e) {
    const std::array<size_t, 4> edofs = {2 * (e + 0) + 0, 2 * (e + 0) + 1,
                                         2 * (e + 1) + 0, 2 * (e + 1) + 1};

    Matrix<real_t, 4, 1> x_elem;
    Matrix<real_t, 4, 1> y_elem;
    Matrix<real_t, 4, 1> z_elem;
    for (size_t a = 0; a < 4; ++a) {
      x_elem(a) = x_prev(edofs[a]);
      y_elem(a) = y_prev(edofs[a]);
      z_elem(a) = z_prev(edofs[a]);
    }

    const Matrix<real_t, 4, 1> bx = coeff * (bending_element_matrix * x_elem);
    const Matrix<real_t, 4, 1> by = coeff * (bending_element_matrix * y_elem);
    const Matrix<real_t, 4, 1> bz = coeff * (bending_element_matrix * z_elem);

    for (size_t a = 0; a < 4; ++a) {
      f_x(edofs[a]) += bx(a);
      f_y(edofs[a]) += by(a);
      f_z(edofs[a]) += bz(a);
    }
  }
}

void EulerBeamInextensibleADDM::apply_dynamic_state_boundary_conditions() {
  if (have_time_dependent_boundary_conditions) {
    apply_time_dependent_boundary_kinematics(time_iter + 1);
    return;
  }

  clear_constrained_dynamic_history();
}

real_t EulerBeamInextensibleADDM::compute_relative_state_update(
    const VectorXd &x_old_iter, const VectorXd &x_new_iter,
    const VectorXd &y_old_iter, const VectorXd &y_new_iter,
    const VectorXd &z_old_iter, const VectorXd &z_new_iter) const {
  const VectorXd dx = x_new_iter - x_old_iter;
  const VectorXd dy = y_new_iter - y_old_iter;
  const VectorXd dz = z_new_iter - z_old_iter;

  const real_t numerator =
      std::sqrt(dx.squaredNorm() + dy.squaredNorm() + dz.squaredNorm());
  const real_t denominator = std::max<real_t>(
      1e-14, std::sqrt(x_old_iter.squaredNorm() + y_old_iter.squaredNorm() +
                       z_old_iter.squaredNorm()));

  return numerator / denominator;
}

real_t EulerBeamInextensibleADDM::compute_max_pq_error() const {
  const real_t res_p = (p - xp).cwiseAbs().maxCoeff();
  const real_t res_q = (q - yp).cwiseAbs().maxCoeff();
  const real_t res_r = (r - zp).cwiseAbs().maxCoeff();
  return std::max(res_p, std::max(res_q, res_r));
}

real_t EulerBeamInextensibleADDM::compute_max_state_update(
    const VectorXd &x_old_iter, const VectorXd &x_new_iter,
    const VectorXd &y_old_iter, const VectorXd &y_new_iter,
    const VectorXd &z_old_iter, const VectorXd &z_new_iter) const {
  const real_t max_dx = (x_new_iter - x_old_iter).cwiseAbs().maxCoeff();
  const real_t max_dy = (y_new_iter - y_old_iter).cwiseAbs().maxCoeff();
  const real_t max_dz = (z_new_iter - z_old_iter).cwiseAbs().maxCoeff();
  return std::max({max_dx, max_dy, max_dz});
}

void EulerBeamInextensibleADDM::update_average_acceleration_state_component(
    const VectorXd &u_old, const VectorXd &u_new, VectorXd &v_hist,
    VectorXd &a_hist, real_t dt) {
  const VectorXd v_old = v_hist;
  v_hist = (2.0 / dt) * (u_new - u_old) - v_old;
  a_hist = (v_hist - v_old) / dt;
}

} // namespace Models
} // namespace ELFF
