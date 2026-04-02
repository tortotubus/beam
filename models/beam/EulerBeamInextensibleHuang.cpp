#include "elff/general/error.hpp"
#include "elff/models/beam/EulerBeamInextensibleHuang.hpp"

#include <cmath>

namespace ELFF {
namespace Models {

EulerBeamInextensibleHuang::EulerBeamInextensibleHuang(real_t length,
                               real_t EI,
                               real_t mu,
                               size_t nodes,
                               EulerBeam::EulerBeamBCs bcs)
  : EulerBeam(length, EI, mu, nodes, bcs)
  , ds(mesh.get_ds())
  , tol_inner(1e-12)
  , last_dt(0.0)
  , initial_velocity_pending(true)
  , implicit_bending(false)
  , X_init(MatX3::Zero(static_cast<Index>(nodes), 3))
  , X_n(MatX3::Zero(static_cast<Index>(nodes), 3))
  , X_nm1(MatX3::Zero(static_cast<Index>(nodes), 3))
  , V_init(MatX3::Zero(static_cast<Index>(nodes), 3))
  , last_T(VectorXd::Zero(static_cast<Index>(nodes > 0 ? nodes - 1 : 0)))
  , bending_matrix(
      MatrixXd::Zero(static_cast<Index>(nodes), static_cast<Index>(nodes)))
{
  dimension = 3;
  validate_boundary_conditions();
  apply_initial_condition();
  bending_matrix = build_bending_matrix();
}

void
EulerBeamInextensibleHuang::set_implicit_bending(bool enabled)
{
  implicit_bending = enabled;
}

bool
EulerBeamInextensibleHuang::get_implicit_bending() const
{
  return implicit_bending;
}

void
EulerBeamInextensibleHuang::solve(std::array<real_t, 3> load)
{
  static_cast<void>(load);
  ELFF_ABORT("EulerBeamInextensibleHuang only implements the dynamic solve(dt, load) "
             "path in this prototype.");
}

void
EulerBeamInextensibleHuang::solve(std::vector<std::array<real_t, 3>> load)
{
  static_cast<void>(load);
  ELFF_ABORT("EulerBeamInextensibleHuang only implements the dynamic solve(dt, load) "
             "path in this prototype.");
}

void
EulerBeamInextensibleHuang::solve(real_t dt, std::array<real_t, 3> load)
{
  ELFF_ASSERT(dt > 0.0, "Time step must be positive.");
  ELFF_ASSERT(mu > 0.0,
              "EulerBeamInextensibleHuang requires mu > 0 for the dynamic solve.");

  if (mesh.get_nodes() < 3) {
    ELFF_ABORT("EulerBeamInextensibleHuang requires at least 3 nodes for the BC1 "
               "finite-difference stencil.");
  }

  if (initial_velocity_pending) {
    X_nm1 = X_n - dt * V_init;
    initial_velocity_pending = false;
  }

  last_dt = dt;

  Vec3 body_force = Vec3::Zero();
  body_force(0) = load[0];
  body_force(1) = load[1];
  body_force(2) = load[2];
  const Vec3 body_accel = body_force / mu;

  const MatX3 X_star = 2.0 * X_n - X_nm1;
  const MatX3 Fb_star = bending_force(X_star);

  // The Huang tension solve is formulated in specific-tension units
  // consistent with acceleration-form forcing. Convert back to a physical
  // tension so the position equation matches the force-balance form used by
  // the other beam solvers.
  last_T = mu * solve_tension(dt, body_accel);
  MatX3 X_np1 = solve_position(last_T, dt, body_force, Fb_star);

  const auto [inext_linf, inext_l2] = compute_inextensibility_error(X_np1);
  ELFF_LOG(time_iter << "\t" << (t + dt) << "\t|g|_inf=" << inext_linf
                     << "\t|g|_l2=" << inext_l2);

  X_nm1 = X_n;
  X_n = X_np1;

  update_mesh();
  time_iter++;
  t += dt;
}

void
EulerBeamInextensibleHuang::solve(real_t dt, std::vector<std::array<real_t, 3>> load)
{
  ELFF_ASSERT(load.size() == mesh.get_nodes(),
              "Nonuniform load must be specified at every node.");

  if (load.empty()) {
    solve(dt, std::array<real_t, 3>{ 0.0, 0.0, 0.0 });
    return;
  }

  const auto& load0 = load.front();
  for (size_t i = 1; i < load.size(); ++i) {
    const auto& li = load[i];
    const real_t diff = std::abs(li[0] - load0[0]) +
                        std::abs(li[1] - load0[1]) +
                        std::abs(li[2] - load0[2]);
    if (diff > tol_inner) {
      ELFF_ABORT("EulerBeamInextensibleHuang prototype currently supports only uniform "
                 "distributed loading.");
    }
  }

  solve(dt, load0);
}

void
EulerBeamInextensibleHuang::apply_initial_condition()
{
  apply_initial_condition(mesh);
}

void
EulerBeamInextensibleHuang::apply_initial_condition(EulerBeamMesh& bmesh)
{
  ELFF_ASSERT(mesh.get_nodes() == bmesh.get_nodes(),
              "Provided mesh must have the same number of nodes.");

  const size_t nodes = mesh.get_nodes();
  auto& centerline = bmesh.get_centerline();
  auto& velocity = bmesh.get_centerline_velocity();

  for (size_t i = 0; i < nodes; ++i) {
    const Index ii = static_cast<Index>(i);
    X_init(ii, 0) = centerline[i][0];
    X_init(ii, 1) = centerline[i][1];
    X_init(ii, 2) = centerline[i][2];
    X_n(ii, 0) = centerline[i][0];
    X_n(ii, 1) = centerline[i][1];
    X_n(ii, 2) = centerline[i][2];
    X_nm1.row(ii) = X_n.row(ii);
    V_init(ii, 0) = velocity[i][0];
    V_init(ii, 1) = velocity[i][1];
    V_init(ii, 2) = velocity[i][2];
  }

  initial_velocity_pending = true;
  update_mesh();
}

const EulerBeamInextensibleHuang::MatX3&
EulerBeamInextensibleHuang::get_initial_positions() const
{
  return X_init;
}

const EulerBeamInextensibleHuang::MatX3&
EulerBeamInextensibleHuang::get_positions() const
{
  return X_n;
}

const VectorXd&
EulerBeamInextensibleHuang::get_last_tension() const
{
  return last_T;
}

Index
EulerBeamInextensibleHuang::bc_index(EulerBeam::EulerBeamBCEnd end) const
{
  if (boundary_conditions.end[0] == end) {
    return 0;
  }
  if (boundary_conditions.end[1] == end) {
    return 1;
  }
  ELFF_ABORT("Boundary condition end mapping is incomplete.");
}

EulerBeam::EulerBeamBCType
EulerBeamInextensibleHuang::bc_type(EulerBeam::EulerBeamBCEnd end) const
{
  return boundary_conditions.type[bc_index(end)];
}

const EulerBeam::EulerBeamBCVals&
EulerBeamInextensibleHuang::bc_vals(EulerBeam::EulerBeamBCEnd end) const
{
  return boundary_conditions.vals[bc_index(end)];
}

bool
EulerBeamInextensibleHuang::is_free(EulerBeam::EulerBeamBCEnd end) const
{
  return bc_type(end) == EulerBeam::free_bc;
}

bool
EulerBeamInextensibleHuang::is_simple(EulerBeam::EulerBeamBCEnd end) const
{
  return bc_type(end) == EulerBeam::simple_bc;
}

bool
EulerBeamInextensibleHuang::is_clamped(EulerBeam::EulerBeamBCEnd end) const
{
  return bc_type(end) == EulerBeam::clamped_bc;
}

bool
EulerBeamInextensibleHuang::is_supported(EulerBeam::EulerBeamBCEnd end) const
{
  return is_simple(end) || is_clamped(end);
}

bool
EulerBeamInextensibleHuang::is_supported_type(EulerBeam::EulerBeamBCType type) const
{
  return type == EulerBeam::free_bc || type == EulerBeam::simple_bc ||
         type == EulerBeam::clamped_bc;
}

void
EulerBeamInextensibleHuang::validate_boundary_conditions() const
{
  ELFF_ASSERT(
    is_supported_type(boundary_conditions.type[0]) &&
      is_supported_type(boundary_conditions.type[1]),
    "EulerBeamInextensibleHuang supports only free_bc, simple_bc, and clamped_bc.");
}

EulerBeamInextensibleHuang::Vec3
EulerBeamInextensibleHuang::endpoint_position(EulerBeam::EulerBeamBCEnd end) const
{
  Vec3 x = Vec3::Zero();
  x(0) = bc_vals(end).position[0];
  x(1) = bc_vals(end).position[1];
  x(2) = bc_vals(end).position[2];
  return x;
}

EulerBeamInextensibleHuang::Vec3
EulerBeamInextensibleHuang::endpoint_tangent(EulerBeam::EulerBeamBCEnd end) const
{
  Vec3 tau = Vec3::Zero();
  tau(0) = bc_vals(end).slope[0];
  tau(1) = bc_vals(end).slope[1];
  tau(2) = bc_vals(end).slope[2];
  ELFF_ASSERT(tau.norm() > tol_inner,
              "Clamped EulerBeamInextensibleHuang boundary requires a nonzero slope.");
  return tau;
}

EulerBeamInextensibleHuang::MatX3
EulerBeamInextensibleHuang::tau_half(const MatX3& X) const
{
  const Index segments = static_cast<Index>(mesh.get_nodes() - 1);
  MatX3 tau = MatX3::Zero(segments, 3);
  for (Index i = 0; i < segments; ++i) {
    tau.row(i) = (X.row(i + 1) - X.row(i)) / ds;
  }
  return tau;
}

std::pair<real_t, real_t>
EulerBeamInextensibleHuang::compute_inextensibility_error(const MatX3& X) const
{
  const MatX3 tau = tau_half(X);
  const Index segments = tau.rows();

  if (segments == 0) {
    return { 0.0, 0.0 };
  }

  real_t max_abs = 0.0;
  real_t sum_sq = 0.0;

  for (Index i = 0; i < segments; ++i) {
    const real_t gi = tau.row(i).squaredNorm() - 1.0;
    const real_t abs_gi = std::abs(gi);
    max_abs = std::max(max_abs, abs_gi);
    sum_sq += gi * gi;
  }

  return { max_abs, std::sqrt(sum_sq / static_cast<real_t>(segments)) };
}

EulerBeamInextensibleHuang::MatX3
EulerBeamInextensibleHuang::dss_nodes(const MatX3& X) const
{
  const Index n = static_cast<Index>(mesh.get_nodes());
  MatX3 dss = MatX3::Zero(n, 3);
  for (Index i = 1; i < n - 1; ++i) {
    dss.row(i) = (X.row(i + 1) - 2.0 * X.row(i) + X.row(i - 1)) / (ds * ds);
  }

  if (is_clamped(EulerBeam::left)) {
    dss.row(0) = (((X.row(1) - X.row(0)) / ds) -
                  endpoint_tangent(EulerBeam::left).transpose()) /
                 (0.5 * ds);
  } else if (!(is_free(EulerBeam::left) || is_simple(EulerBeam::left))) {
    ELFF_ABORT("Unsupported left boundary condition for EulerBeamInextensibleHuang.");
  }

  if (is_clamped(EulerBeam::right)) {
    dss.row(n - 1) = (endpoint_tangent(EulerBeam::right).transpose() -
                      ((X.row(n - 1) - X.row(n - 2)) / ds)) /
                     (0.5 * ds);
  } else if (!(is_free(EulerBeam::right) || is_simple(EulerBeam::right))) {
    ELFF_ABORT("Unsupported right boundary condition for EulerBeamInextensibleHuang.");
  }

  return dss;
}

EulerBeamInextensibleHuang::MatX3
EulerBeamInextensibleHuang::bending_force(const MatX3& X) const
{
  const Index n = static_cast<Index>(mesh.get_nodes());
  MatX3 dss = dss_nodes(X);
  MatX3 Fb = MatX3::Zero(n, 3);

  if (is_free(EulerBeam::left) && n >= 3) {
    Fb.row(0) = -EI * (dss.row(2) - dss.row(1)) / (ds * ds);
  }

  for (Index i = 1; i < n - 1; ++i) {
    Fb.row(i) =
      -EI * (dss.row(i + 1) - 2.0 * dss.row(i) + dss.row(i - 1)) / (ds * ds);
  }

  if (is_free(EulerBeam::right) && n >= 3) {
    Fb.row(n - 1) = -EI * (dss.row(n - 3) - dss.row(n - 2)) / (ds * ds);
  }

  return Fb;
}

MatrixXd
EulerBeamInextensibleHuang::build_bending_matrix() const
{
  const Index n = static_cast<Index>(mesh.get_nodes());
  MatrixXd B = MatrixXd::Zero(n, n);

  for (Index j = 0; j < n; ++j) {
    MatX3 X_basis = MatX3::Zero(n, 3);
    X_basis(j, 0) = 1.0;
    B.col(j) = bending_force(X_basis).col(0);
  }

  return B;
}

VectorXd
EulerBeamInextensibleHuang::solve_tension(real_t dt, const Vec3& body_force) const
{
  const Index segments = static_cast<Index>(mesh.get_nodes() - 1);
  const Index last = segments - 1;
  ELFF_ASSERT(mu > 0.0,
              "EulerBeamInextensibleHuang requires mu > 0 for the tension solve.");

  MatX3 X_star = 2.0 * X_n - X_nm1;
  MatX3 U_n = (X_n - X_nm1) / dt;

  MatX3 tau_star = tau_half(X_star);
  MatX3 tau_n = tau_half(X_n);
  MatX3 tau_nm1 = tau_half(X_nm1);
  MatX3 Fb_star = bending_force(X_star) / mu;

  VectorXd corr = VectorXd::Zero(segments);
  VectorXd velsq = VectorXd::Zero(segments);
  for (Index i = 0; i < segments; ++i) {
    corr(i) = (1.0 - 2.0 * tau_n.row(i).squaredNorm() +
               tau_nm1.row(i).squaredNorm()) /
              (2.0 * dt * dt);
  }

  MatX3 Du_n = tau_half(U_n);
  for (Index i = 0; i < segments; ++i) {
    velsq(i) = Du_n.row(i).squaredNorm();
  }

  MatrixXd A = MatrixXd::Zero(segments, segments);
  VectorXd b = VectorXd::Zero(segments);

  for (Index i = 1; i < segments - 1; ++i) {
    A(i, i - 1) = tau_star.row(i).dot(tau_star.row(i - 1)) / (ds * ds);
    A(i, i) = -2.0 * tau_star.row(i).squaredNorm() / (ds * ds);
    A(i, i + 1) = tau_star.row(i).dot(tau_star.row(i + 1)) / (ds * ds);
    b(i) =
      corr(i) - velsq(i) -
      tau_star.row(i).dot((Fb_star.row(i + 1) - Fb_star.row(i)).transpose()) /
        ds;
  }

  if (is_free(EulerBeam::left)) {
    A(0, 0) = -3.0 * tau_star.row(0).squaredNorm() / (ds * ds);
    if (segments > 1) {
      A(0, 1) = tau_star.row(0).dot(tau_star.row(1)) / (ds * ds);
    }
    b(0) =
      corr(0) - velsq(0) -
      tau_star.row(0).dot((Fb_star.row(1) - Fb_star.row(0)).transpose()) / ds;
  } else {
    A(0, 0) = -tau_star.row(0).squaredNorm() / (ds * ds);
    if (segments > 1) {
      A(0, 1) = tau_star.row(0).dot(tau_star.row(1)) / (ds * ds);
    }
    b(0) = corr(0) - velsq(0) +
           tau_star.row(0).dot((-Fb_star.row(1)).transpose()) / ds -
           tau_star.row(0).dot(body_force.transpose()) / ds;
  }

  if (is_free(EulerBeam::right)) {
    if (segments > 1) {
      A(last, last - 1) =
        tau_star.row(last).dot(tau_star.row(last - 1)) / (ds * ds);
    }
    A(last, last) = -3.0 * tau_star.row(last).squaredNorm() / (ds * ds);
    b(last) =
      corr(last) - velsq(last) +
      tau_star.row(last).dot(
        (Fb_star.row(segments - 1) - Fb_star.row(segments)).transpose()) /
        ds;
  } else {
    if (segments > 1) {
      A(last, last - 1) =
        tau_star.row(last).dot(tau_star.row(last - 1)) / (ds * ds);
    }
    A(last, last) = -tau_star.row(last).squaredNorm() / (ds * ds);
    b(last) = corr(last) - velsq(last) +
              tau_star.row(last).dot(Fb_star.row(last).transpose()) / ds +
              tau_star.row(last).dot(body_force.transpose()) / ds;
  }

  return A.partialPivLu().solve(b);
}

EulerBeamInextensibleHuang::MatX3
EulerBeamInextensibleHuang::solve_position(const VectorXd& T_half,
                               real_t dt,
                               const Vec3& body_force,
                               const MatX3& Fb_star) const
{
  const Index n = static_cast<Index>(mesh.get_nodes());
  ELFF_ASSERT(mu > 0.0,
              "EulerBeamInextensibleHuang requires mu > 0 for the position solve.");

  MatX3 rhs = mu * (2.0 * X_n - X_nm1) / (dt * dt);
  for (Index i = 0; i < n; ++i) {
    rhs.row(i) += body_force.transpose();
    if (!implicit_bending) {
      rhs.row(i) += Fb_star.row(i);
    }
  }

  MatX3 X_np1 = MatX3::Zero(n, 3);
  const Vec3 x_left = endpoint_position(EulerBeam::left);
  const Vec3 x_right = endpoint_position(EulerBeam::right);
  const Index last_segment = static_cast<Index>(T_half.size() - 1);

  for (Index comp = 0; comp < 3; ++comp) {
    MatrixXd A = MatrixXd::Zero(n, n);
    VectorXd b = rhs.col(comp);

    if (is_free(EulerBeam::left)) {
      A(0, 0) = mu / (dt * dt) + 2.0 * T_half(0) / (ds * ds);
      A(0, 1) = -2.0 * T_half(0) / (ds * ds);
    }

    for (Index i = 1; i < n - 1; ++i) {
      A(i, i - 1) = -T_half(i - 1) / (ds * ds);
      A(i, i) = mu / (dt * dt) + (T_half(i - 1) + T_half(i)) / (ds * ds);
      A(i, i + 1) = -T_half(i) / (ds * ds);
    }

    if (is_free(EulerBeam::right)) {
      A(n - 1, n - 2) = -2.0 * T_half(last_segment) / (ds * ds);
      A(n - 1, n - 1) =
        mu / (dt * dt) + 2.0 * T_half(last_segment) / (ds * ds);
    }

    if (implicit_bending) {
      A -= bending_matrix;
    }

    if (is_supported(EulerBeam::left)) {
      A.row(0).setZero();
      A(0, 0) = 1.0;
      b(0) = x_left(comp);
    }

    if (is_supported(EulerBeam::right)) {
      A.row(n - 1).setZero();
      A(n - 1, n - 1) = 1.0;
      b(n - 1) = x_right(comp);
    }

    X_np1.col(comp) = A.partialPivLu().solve(b);
  }

  return X_np1;
}

void
EulerBeamInextensibleHuang::update_mesh()
{
  const size_t nodes = mesh.get_nodes();
  auto& centerline = mesh.get_centerline();
  auto& velocity = mesh.get_centerline_velocity();
  auto& slope = mesh.get_slope();

  for (size_t i = 0; i < nodes; ++i) {
    const Index ii = static_cast<Index>(i);

    centerline[i][0] = X_n(ii, 0);
    centerline[i][1] = X_n(ii, 1);
    centerline[i][2] = X_n(ii, 2);

    if (last_dt > 0.0) {
      velocity[i][0] = (X_n(ii, 0) - X_nm1(ii, 0)) / last_dt;
      velocity[i][1] = (X_n(ii, 1) - X_nm1(ii, 1)) / last_dt;
      velocity[i][2] = (X_n(ii, 2) - X_nm1(ii, 2)) / last_dt;
    } else {
      velocity[i][0] = V_init(ii, 0);
      velocity[i][1] = V_init(ii, 1);
      velocity[i][2] = V_init(ii, 2);
    }

    if (i == 0) {
      slope[i][0] = (X_n(1, 0) - X_n(0, 0)) / ds;
      slope[i][1] = (X_n(1, 1) - X_n(0, 1)) / ds;
      slope[i][2] = (X_n(1, 2) - X_n(0, 2)) / ds;
    } else if (i + 1 == nodes) {
      slope[i][0] = (X_n(ii, 0) - X_n(ii - 1, 0)) / ds;
      slope[i][1] = (X_n(ii, 1) - X_n(ii - 1, 1)) / ds;
      slope[i][2] = (X_n(ii, 2) - X_n(ii - 1, 2)) / ds;
    } else {
      slope[i][0] = (X_n(ii + 1, 0) - X_n(ii - 1, 0)) / (2.0 * ds);
      slope[i][1] = (X_n(ii + 1, 1) - X_n(ii - 1, 1)) / (2.0 * ds);
      slope[i][2] = (X_n(ii + 1, 2) - X_n(ii - 1, 2)) / (2.0 * ds);
    }
  }
}

} // namespace Models
} // namespace ELFF
