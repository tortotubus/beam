// ============================================================
// EulerBeamInextensibleIndex1.cpp
// ============================================================
#include "elff/models/beam/EulerBeamInextensibleIndex1.hpp"

namespace ELFF {
namespace Models {

namespace {
constexpr std::array<real_t, 3> xi_q = { 0.1127016654, 0.5, 0.8872983346 };
constexpr std::array<real_t, 3>  w_q = { 0.2777777778, 0.4444444444, 0.2777777778 };
} // namespace

// -----------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------
EulerBeamInextensibleIndex1::EulerBeamInextensibleIndex1(
  real_t                  length,
  real_t                  EI,
  real_t                  mu,
  size_t                  n_nodes,
  EulerBeam::EulerBeamBCs bcs)
  : EulerBeam(length, EI, mu, n_nodes, bcs)
  , elements(n_nodes - 1)
  , nodes(n_nodes)
  , ds(mesh.get_ds())
  , ndof_x(2 * n_nodes)
  , ndof_y(2 * n_nodes)
  , ndof_z(2 * n_nodes)
  , ndof_l(n_nodes)
  , offset_x(0)
  , offset_y(ndof_x)
  , offset_z(ndof_x + ndof_y)
  , offset_l(ndof_x + ndof_y + ndof_z)
  , ndof(ndof_x + ndof_y + ndof_z)
  , u(VectorXd::Zero(ndof))
  , lambda(VectorXd::Zero(ndof_l))
  , u_prev(VectorXd::Zero(ndof))
  , v_prev(VectorXd::Zero(ndof))
  , a_prev(VectorXd::Zero(ndof))
  , K_elastic(ndof, ndof)
  , saddle_mat(ndof + ndof_l, ndof + ndof_l)
  , saddle_rhs(VectorXd::Zero(ndof + ndof_l))
  , tol_position_drift(1e-8)
  , tol_velocity_drift(1e-8)
  , max_projection_iter(20)
{
  apply_initial_condition(mesh);
  assemble_elastic_stiffness();
  u_prev = u;
}

// -----------------------------------------------------------------------
// Initial conditions
// -----------------------------------------------------------------------
void EulerBeamInextensibleIndex1::solve()
{
  ELFF_ABORT("EulerBeamInextensibleIndex1 does not implement a static solve.\n");
}

void EulerBeamInextensibleIndex1::solve(std::array<real_t, 3> load)
{
  static_cast<void>(load);
  ELFF_ABORT("EulerBeamInextensibleIndex1 does not implement a static solve.\n");
}

void EulerBeamInextensibleIndex1::solve(
  std::vector<std::array<real_t, 3>> load)
{
  static_cast<void>(load);
  ELFF_ABORT("EulerBeamInextensibleIndex1 does not implement a static solve.\n");
}

void EulerBeamInextensibleIndex1::apply_initial_condition()
{
  apply_initial_condition(mesh);
}

void EulerBeamInextensibleIndex1::apply_initial_condition(
  EulerBeamMesh& bmesh)
{
  ELFF_ASSERT(nodes == bmesh.get_nodes(),
              "Provided mesh must have same node count as the previous mesh.\n");

  const auto centerline = bmesh.get_centerline();
  const auto slopes = bmesh.get_slope();
  const auto velocity = bmesh.get_centerline_velocity();

  for (size_t ni = 0; ni < nodes; ++ni) {
    u(offset_x + 2 * ni + 0) = centerline[ni][0];
    u(offset_x + 2 * ni + 1) = slopes[ni][0];
    u(offset_y + 2 * ni + 0) = centerline[ni][1];
    u(offset_y + 2 * ni + 1) = slopes[ni][1];
    u(offset_z + 2 * ni + 0) = centerline[ni][2];
    u(offset_z + 2 * ni + 1) = slopes[ni][2];

    v_prev(offset_x + 2 * ni) = velocity[ni][0];
    v_prev(offset_y + 2 * ni) = velocity[ni][1];
    v_prev(offset_z + 2 * ni) = velocity[ni][2];
  }

  lambda.setZero();
  u_prev = u;
  a_prev.setZero();
  time_iter = 0;
  t = 0.0;
  update_mesh();
}

// -----------------------------------------------------------------------
// Public solve entry points
// -----------------------------------------------------------------------
void EulerBeamInextensibleIndex1::solve(real_t              dt,
                                                std::array<real_t, 3> load)
{
  if (!(dt > 0.0))
    throw std::runtime_error("Index-1 Newmark: dt must be > 0");

  constexpr real_t alpha = 0.0;
  constexpr real_t gamma = 0.5 - alpha;
  const real_t     beta  = 0.25 * (1.0 - alpha) * (1.0 - alpha);

  // Assemble constraint Jacobian B and velocity-level rhs g
  // from the current configuration u and velocity v_prev.
  SparseMatrix<real_t> B;
  VectorXd             g(ndof_l);
  assemble_B_and_g(B, g);

  // External force vector
  const VectorXd f_ext = assemble_f_ext(load);

  // Build and solve the saddle-point system
  assemble_saddle_system(dt, beta, f_ext, B, g);
  apply_saddle_boundary_conditions();

  SparseLU<SparseMatrix<real_t>, COLAMDOrdering<int>> solver;
  solver.compute(saddle_mat);
  if (solver.info() != Success)
    ELFF_ABORT("EulerBeamInextensibleIndex1::solve: LU factorization failed\n");

  const VectorXd sol = solver.solve(saddle_rhs);
  if (solver.info() != Success)
    ELFF_ABORT("EulerBeamInextensibleIndex1::solve: back-solve failed\n");

  // Extract displacement and physical tension
  u      = sol.head(ndof);
  lambda = sol.tail(ndof_l);

  // Newmark velocity / acceleration update
  update_newmark_state(dt, beta, gamma);
  apply_dynamic_state_boundary_conditions();

  // Drift correction: position first, then velocity with updated geometry
  const real_t pos_drift = assemble_constraint_residual_vector().norm();
  if (pos_drift > tol_position_drift) {
    project_position_onto_constraint();
    // Reassemble B at corrected position for velocity projection
    SparseMatrix<real_t> B_corr;
    VectorXd             g_dummy(ndof_l);
    assemble_B_and_g(B_corr, g_dummy);
    project_velocity_onto_constraint(B_corr);
  } else if ((B * v_prev).norm() > tol_velocity_drift) {
    project_velocity_onto_constraint(B);
  }

  u_prev = u;

  update_mesh();
  ++time_iter;
  t += dt;
}

void EulerBeamInextensibleIndex1::solve(
  real_t                              dt,
  std::vector<std::array<real_t, 3>> load)
{
  ELFF_ASSERT(load.size() == nodes,
              "Load vector size must equal node count.\n");

  if (!(dt > 0.0))
    throw std::runtime_error("Index-1 Newmark: dt must be > 0");

  constexpr real_t alpha = 0.0;
  constexpr real_t gamma = 0.5 - alpha;
  const real_t     beta  = 0.25 * (1.0 - alpha) * (1.0 - alpha);

  SparseMatrix<real_t> B;
  VectorXd             g(ndof_l);
  assemble_B_and_g(B, g);

  const VectorXd f_ext = assemble_f_ext(load);

  assemble_saddle_system(dt, beta, f_ext, B, g);
  apply_saddle_boundary_conditions();

  SparseLU<SparseMatrix<real_t>, COLAMDOrdering<int>> solver;
  solver.compute(saddle_mat);
  if (solver.info() != Success)
    ELFF_ABORT("EulerBeamInextensibleIndex1::solve: LU factorization failed\n");

  const VectorXd sol = solver.solve(saddle_rhs);
  if (solver.info() != Success)
    ELFF_ABORT("EulerBeamInextensibleIndex1::solve: back-solve failed\n");

  u      = sol.head(ndof);
  lambda = sol.tail(ndof_l);

  update_newmark_state(dt, beta, gamma);
  apply_dynamic_state_boundary_conditions();

  const real_t pos_drift = assemble_constraint_residual_vector().norm();
  if (pos_drift > tol_position_drift) {
    project_position_onto_constraint();
    SparseMatrix<real_t> B_corr;
    VectorXd             g_dummy(ndof_l);
    assemble_B_and_g(B_corr, g_dummy);
    project_velocity_onto_constraint(B_corr);
  } else if ((B * v_prev).norm() > tol_velocity_drift) {
    project_velocity_onto_constraint(B);
  }

  u_prev = u;

  update_mesh();
  ++time_iter;
  t += dt;
}

// -----------------------------------------------------------------------
// Elastic stiffness (constant, assembled once)
//   K_elastic = EI * int H_a'' H_b'' ds   (same for x, y, z - decoupled)
// -----------------------------------------------------------------------
void EulerBeamInextensibleIndex1::assemble_elastic_stiffness()
{
  using Tpl = Triplet<real_t>;
  std::vector<Tpl> triplets;
  triplets.reserve(elements * 3 * 16);

  for (size_t e = 0; e < elements; ++e) {
    const size_t ix[4] = { offset_x + 2 * e,       offset_x + 2 * e + 1,
                            offset_x + 2 * (e + 1), offset_x + 2 * (e + 1) + 1 };
    const size_t iy[4] = { offset_y + 2 * e,       offset_y + 2 * e + 1,
                            offset_y + 2 * (e + 1), offset_y + 2 * (e + 1) + 1 };
    const size_t iz[4] = { offset_z + 2 * e,       offset_z + 2 * e + 1,
                            offset_z + 2 * (e + 1), offset_z + 2 * (e + 1) + 1 };

    real_t K4e[4][4] = {};

    for (size_t qi = 0; qi < 3; ++qi) {
      const auto   ddH = ELFF::FEM::CubicHermite<real_t>::second_derivs(xi_q[qi], ds);
      const real_t w   = w_q[qi];
      for (size_t a = 0; a < 4; ++a)
        for (size_t b = 0; b < 4; ++b)
          K4e[a][b] += ddH[a] * ddH[b] * w * ds;
    }

    for (size_t a = 0; a < 4; ++a) {
      for (size_t b = 0; b < 4; ++b) {
        const real_t v = EI * K4e[a][b];
        triplets.emplace_back(ix[a], ix[b], v);
        triplets.emplace_back(iy[a], iy[b], v);
        triplets.emplace_back(iz[a], iz[b], v);
      }
    }
  }

  K_elastic.resize(ndof, ndof);
  K_elastic.setFromTriplets(triplets.begin(), triplets.end());
  K_elastic.makeCompressed();
}

// -----------------------------------------------------------------------
// Constraint Jacobian B and velocity-level RHS g
//
//   B_{a,i_c} = int eta_a  *  r'_c  *  dH_i  ds     c in {x,y,z}
//   g_a       = -int eta_a * ||r_dot'||^2 ds
//
// B has size ndof_l x ndof.  g has size ndof_l.
// -----------------------------------------------------------------------
void EulerBeamInextensibleIndex1::assemble_B_and_g(
  SparseMatrix<real_t>& B,
  VectorXd&             g) const
{
  using Tpl = Triplet<real_t>;
  std::vector<Tpl> triplets;
  triplets.reserve(elements * 2 * 12);
  g = VectorXd::Zero(ndof_l);

  for (size_t e = 0; e < elements; ++e) {
    const size_t n0 = e, n1 = e + 1;
    const size_t ix[4] = { offset_x + 2 * n0, offset_x + 2 * n0 + 1,
                            offset_x + 2 * n1, offset_x + 2 * n1 + 1 };
    const size_t iy[4] = { offset_y + 2 * n0, offset_y + 2 * n0 + 1,
                            offset_y + 2 * n1, offset_y + 2 * n1 + 1 };
    const size_t iz[4] = { offset_z + 2 * n0, offset_z + 2 * n0 + 1,
                            offset_z + 2 * n1, offset_z + 2 * n1 + 1 };
    const size_t il[2] = { n0, n1 };

    for (size_t qi = 0; qi < 3; ++qi) {
      const real_t xi_q_val = xi_q[qi];
      const real_t w        = w_q[qi];

      const auto dH = ELFF::FEM::CubicHermite<real_t>::derivs(xi_q_val, ds);
      const auto M  = ELFF::FEM::LinearShape<real_t>::values(xi_q_val);

      // Current slope r' and slope velocity r_dot' at quadrature point
      real_t xp = 0., yp = 0., zp = 0.;
      real_t vxp = 0., vyp = 0., vzp = 0.;
      for (size_t i = 0; i < 4; ++i) {
        xp  += dH[i] * u(ix[i]);
        yp  += dH[i] * u(iy[i]);
        zp  += dH[i] * u(iz[i]);
        vxp += dH[i] * v_prev(ix[i]);
        vyp += dH[i] * v_prev(iy[i]);
        vzp += dH[i] * v_prev(iz[i]);
      }

      // Velocity-level RHS: -||r_dot'||^2
      const real_t g_q = -(vxp * vxp + vyp * vyp + vzp * vzp);

      for (size_t a = 0; a < 2; ++a) {
        const real_t Ma_w_ds = M[a] * w * ds;

        // g contribution
        g(il[a]) += g_q * Ma_w_ds;

        // B entries: B[l_a][disp_i] = int M_a * r'_c * dH_i ds
        for (size_t i = 0; i < 4; ++i) {
          const real_t c = Ma_w_ds * dH[i];
          triplets.emplace_back(il[a], ix[i], xp * c);
          triplets.emplace_back(il[a], iy[i], yp * c);
          triplets.emplace_back(il[a], iz[i], zp * c);
        }
      }
    }
  }

  B.resize(ndof_l, ndof);
  B.setFromTriplets(triplets.begin(), triplets.end());
  B.makeCompressed();
}

// -----------------------------------------------------------------------
// Position constraint residual vector (used for drift detection)
//   c_a = int eta_a * (||r'||^2 - 1) ds
// -----------------------------------------------------------------------
VectorXd EulerBeamInextensibleIndex1::assemble_constraint_residual_vector()
  const
{
  VectorXd c = VectorXd::Zero(ndof_l);

  for (size_t e = 0; e < elements; ++e) {
    const size_t n0 = e, n1 = e + 1;
    const size_t ix[4] = { offset_x + 2 * n0, offset_x + 2 * n0 + 1,
                            offset_x + 2 * n1, offset_x + 2 * n1 + 1 };
    const size_t iy[4] = { offset_y + 2 * n0, offset_y + 2 * n0 + 1,
                            offset_y + 2 * n1, offset_y + 2 * n1 + 1 };
    const size_t iz[4] = { offset_z + 2 * n0, offset_z + 2 * n0 + 1,
                            offset_z + 2 * n1, offset_z + 2 * n1 + 1 };

    for (size_t qi = 0; qi < 3; ++qi) {
      const real_t xi_q_val = xi_q[qi];
      const real_t w        = w_q[qi];

      const auto dH = ELFF::FEM::CubicHermite<real_t>::derivs(xi_q_val, ds);
      const auto M  = ELFF::FEM::LinearShape<real_t>::values(xi_q_val);

      real_t xp = 0., yp = 0., zp = 0.;
      for (size_t i = 0; i < 4; ++i) {
        xp += dH[i] * u(ix[i]);
        yp += dH[i] * u(iy[i]);
        zp += dH[i] * u(iz[i]);
      }

      const real_t S = xp * xp + yp * yp + zp * zp - 1.0;

      for (size_t a = 0; a < 2; ++a)
        c(n0 + a) += S * M[a] * w * ds;
    }
  }

  return c;
}

// -----------------------------------------------------------------------
// External force assembly - uniform body load
// -----------------------------------------------------------------------
VectorXd EulerBeamInextensibleIndex1::assemble_f_ext(
  std::array<real_t, 3> load) const
{
  VectorXd f = VectorXd::Zero(ndof);

  for (size_t e = 0; e < elements; ++e) {
    const size_t ix[4] = { offset_x + 2 * e,       offset_x + 2 * e + 1,
                            offset_x + 2 * (e + 1), offset_x + 2 * (e + 1) + 1 };
    const size_t iy[4] = { offset_y + 2 * e,       offset_y + 2 * e + 1,
                            offset_y + 2 * (e + 1), offset_y + 2 * (e + 1) + 1 };
    const size_t iz[4] = { offset_z + 2 * e,       offset_z + 2 * e + 1,
                            offset_z + 2 * (e + 1), offset_z + 2 * (e + 1) + 1 };

    for (size_t qi = 0; qi < 3; ++qi) {
      const auto   H = ELFF::FEM::CubicHermite<real_t>::values(xi_q[qi], ds);
      const real_t w = w_q[qi];
      for (size_t a = 0; a < 4; ++a) {
        f(ix[a]) += load[0] * H[a] * w * ds;
        f(iy[a]) += load[1] * H[a] * w * ds;
        f(iz[a]) += load[2] * H[a] * w * ds;
      }
    }
  }

  // Point forces / torques from boundary conditions
  for (size_t bi = 0; bi < 2; ++bi) {
    const size_t          ni     = (boundary_conditions.end[bi] == left) ? 0 : nodes - 1;
    const EulerBeamBCType bctype = boundary_conditions.type[bi];
    const EulerBeamBCVals bcvals = boundary_conditions.vals[bi];

    if (bctype == point_force_bc) {
      f(offset_x + 2 * ni) += bcvals.force[0];
      f(offset_y + 2 * ni) += bcvals.force[1];
      f(offset_z + 2 * ni) += bcvals.force[2];
    } else if (bctype == point_torque_bc) {
      f(offset_x + 2 * ni + 1) += bcvals.torque[0];
      f(offset_y + 2 * ni + 1) += bcvals.torque[1];
      f(offset_z + 2 * ni + 1) += bcvals.torque[2];
    }
  }

  return f;
}

// -----------------------------------------------------------------------
// External force assembly - nodal load vector
// -----------------------------------------------------------------------
VectorXd EulerBeamInextensibleIndex1::assemble_f_ext(
  const std::vector<std::array<real_t, 3>>& load) const
{
  VectorXd f = VectorXd::Zero(ndof);

  for (size_t e = 0; e < elements; ++e) {
    const size_t ix[4] = { offset_x + 2 * e,       offset_x + 2 * e + 1,
                            offset_x + 2 * (e + 1), offset_x + 2 * (e + 1) + 1 };
    const size_t iy[4] = { offset_y + 2 * e,       offset_y + 2 * e + 1,
                            offset_y + 2 * (e + 1), offset_y + 2 * (e + 1) + 1 };
    const size_t iz[4] = { offset_z + 2 * e,       offset_z + 2 * e + 1,
                            offset_z + 2 * (e + 1), offset_z + 2 * (e + 1) + 1 };

    for (size_t qi = 0; qi < 3; ++qi) {
      const auto   H = ELFF::FEM::CubicHermite<real_t>::values(xi_q[qi], ds);
      const auto   M = ELFF::FEM::LinearShape<real_t>::values(xi_q[qi]);
      const real_t w = w_q[qi];

      // Linearly interpolate nodal loads across element
      const real_t fx = M[0] * load[e][0] + M[1] * load[e + 1][0];
      const real_t fy = M[0] * load[e][1] + M[1] * load[e + 1][1];
      const real_t fz = M[0] * load[e][2] + M[1] * load[e + 1][2];

      for (size_t a = 0; a < 4; ++a) {
        f(ix[a]) += fx * H[a] * w * ds;
        f(iy[a]) += fy * H[a] * w * ds;
        f(iz[a]) += fz * H[a] * w * ds;
      }
    }
  }

  for (size_t bi = 0; bi < 2; ++bi) {
    const size_t          ni     = (boundary_conditions.end[bi] == left) ? 0 : nodes - 1;
    const EulerBeamBCType bctype = boundary_conditions.type[bi];
    const EulerBeamBCVals bcvals = boundary_conditions.vals[bi];

    if (bctype == point_force_bc) {
      f(offset_x + 2 * ni) += bcvals.force[0];
      f(offset_y + 2 * ni) += bcvals.force[1];
      f(offset_z + 2 * ni) += bcvals.force[2];
    } else if (bctype == point_torque_bc) {
      f(offset_x + 2 * ni + 1) += bcvals.torque[0];
      f(offset_y + 2 * ni + 1) += bcvals.torque[1];
      f(offset_z + 2 * ni + 1) += bcvals.torque[2];
    }
  }

  return f;
}

// -----------------------------------------------------------------------
// Saddle-point system assembly
//
// [ K_eff   B^T ] [ u^{n+1} ]   [ f_ext + coeff*M*u_tilde ]
// [ B       0   ] [ lambda  ] = [ beta*dt^2*g + B*u_tilde  ]
//
// K_eff = K_elastic + coeff * M_lumped  (coeff = 1/beta/dt^2)
// u_tilde = u^n + dt*v^n + dt^2*(0.5-beta)*a^n
// -----------------------------------------------------------------------
void EulerBeamInextensibleIndex1::assemble_saddle_system(
  real_t                      dt,
  real_t                      beta,
  const VectorXd&             f_ext,
  const SparseMatrix<real_t>& B,
  const VectorXd&             g)
{
  const size_t total_dof = ndof + ndof_l;
  const real_t coeff     = 1.0 / (beta * dt * dt);

  const VectorXd u_tilde =
    u_prev + dt * v_prev + dt * dt * (0.5 - beta) * a_prev;

  using Tpl = Triplet<real_t>;
  std::vector<Tpl> triplets;
  triplets.reserve(K_elastic.nonZeros() + 2 * B.nonZeros() + 3 * nodes);

  // Top-left block: K_elastic + coeff * M_lumped (diagonal)
  for (int k = 0; k < K_elastic.outerSize(); ++k)
    for (SparseMatrix<real_t>::InnerIterator it(K_elastic, k); it; ++it)
      triplets.emplace_back((int)it.row(), (int)it.col(), it.value());

  for (size_t n = 0; n < nodes; ++n) {
    const real_t w_node = (n == 0 || n == nodes - 1) ? 0.5 * ds : ds;
    const real_t m      = mu * w_node * coeff;
    triplets.emplace_back(offset_x + 2 * n, offset_x + 2 * n, m);
    triplets.emplace_back(offset_y + 2 * n, offset_y + 2 * n, m);
    triplets.emplace_back(offset_z + 2 * n, offset_z + 2 * n, m);
  }

  // Top-right block: B^T
  for (int k = 0; k < B.outerSize(); ++k)
    for (SparseMatrix<real_t>::InnerIterator it(B, k); it; ++it)
      triplets.emplace_back((int)it.col(), (int)(ndof + it.row()), it.value());

  // Bottom-left block: B
  for (int k = 0; k < B.outerSize(); ++k)
    for (SparseMatrix<real_t>::InnerIterator it(B, k); it; ++it)
      triplets.emplace_back((int)(ndof + it.row()), (int)it.col(), it.value());

  saddle_mat.resize(total_dof, total_dof);
  saddle_mat.setFromTriplets(triplets.begin(), triplets.end());

  // RHS
  saddle_rhs.resize(total_dof);
  saddle_rhs.setZero();

  // Top: f_ext + coeff * M_lumped * u_tilde
  saddle_rhs.head(ndof) = f_ext;
  for (size_t n = 0; n < nodes; ++n) {
    const real_t w_node  = (n == 0 || n == nodes - 1) ? 0.5 * ds : ds;
    const real_t m_coeff = mu * w_node * coeff;
    const size_t ix = offset_x + 2 * n;
    const size_t iy = offset_y + 2 * n;
    const size_t iz = offset_z + 2 * n;
    saddle_rhs(ix) += m_coeff * u_tilde(ix);
    saddle_rhs(iy) += m_coeff * u_tilde(iy);
    saddle_rhs(iz) += m_coeff * u_tilde(iz);
  }

  // Bottom: beta*dt^2*g + B*u_tilde
  saddle_rhs.tail(ndof_l) = beta * dt * dt * g + B * u_tilde;

  saddle_mat.makeCompressed();
}

// -----------------------------------------------------------------------
// Boundary conditions on the saddle-point system.
//
// Displacement BCs (simple, clamped): eliminate column d by subtracting
// its contribution from the RHS, zero row d, set diagonal to 1.
//
// Lambda at free ends: zero the lambda row, set diagonal to 1, RHS = 0.
// At free ends the physical tension is zero (natural BC for tension).
// -----------------------------------------------------------------------
void EulerBeamInextensibleIndex1::apply_saddle_boundary_conditions()
{
  // Collect displacement-constrained DOFs and their prescribed values
  std::vector<std::pair<size_t, real_t>> disp_bcs;

  for (size_t bi = 0; bi < 2; ++bi) {
    const size_t          ni     = (boundary_conditions.end[bi] == left) ? 0 : nodes - 1;
    const EulerBeamBCType bctype = boundary_conditions.type[bi];
    const EulerBeamBCVals bcvals = boundary_conditions.vals[bi];

    switch (bctype) {
      case simple_bc:
        disp_bcs.push_back({ offset_x + 2 * ni, bcvals.position[0] });
        disp_bcs.push_back({ offset_y + 2 * ni, bcvals.position[1] });
        disp_bcs.push_back({ offset_z + 2 * ni, bcvals.position[2] });
        break;
      case clamped_bc:
        disp_bcs.push_back({ offset_x + 2 * ni,     bcvals.position[0] });
        disp_bcs.push_back({ offset_x + 2 * ni + 1, bcvals.slope[0] });
        disp_bcs.push_back({ offset_y + 2 * ni,     bcvals.position[1] });
        disp_bcs.push_back({ offset_y + 2 * ni + 1, bcvals.slope[1] });
        disp_bcs.push_back({ offset_z + 2 * ni,     bcvals.position[2] });
        disp_bcs.push_back({ offset_z + 2 * ni + 1, bcvals.slope[2] });
        break;
      default:
        break;
    }
  }

  for (auto& [d, val] : disp_bcs) {
    // Column d: subtract contribution from RHS then zero (CSC-efficient)
    for (SparseMatrix<real_t>::InnerIterator it(saddle_mat, d); it; ++it) {
      if ((size_t)it.row() != d) {
        saddle_rhs(it.row()) -= it.value() * val;
        it.valueRef() = 0.0;
      }
    }
    // Row d: zero all entries (requires full scan in CSC format)
    for (int col = 0; col < saddle_mat.outerSize(); ++col) {
      for (SparseMatrix<real_t>::InnerIterator it(saddle_mat, col); it; ++it) {
        if ((size_t)it.row() == d && (size_t)it.col() != d)
          it.valueRef() = 0.0;
      }
    }
    saddle_mat.coeffRef(d, d) = 1.0;
    saddle_rhs(d)             = val;
  }

  // Lambda = 0 at free ends (zero tension, natural BC)
  for (size_t bi = 0; bi < 2; ++bi) {
    if (boundary_conditions.type[bi] != free_bc) continue;
    const size_t ni    = (boundary_conditions.end[bi] == left) ? 0 : nodes - 1;
    const size_t l_dof = ndof + ni;

    for (int col = 0; col < saddle_mat.outerSize(); ++col) {
      for (SparseMatrix<real_t>::InnerIterator it(saddle_mat, col); it; ++it) {
        if ((size_t)it.row() == l_dof && (size_t)it.col() != l_dof)
          it.valueRef() = 0.0;
      }
    }
    saddle_mat.coeffRef(l_dof, l_dof) = 1.0;
    saddle_rhs(l_dof)                 = 0.0;
  }

  saddle_mat.makeCompressed();
}

// -----------------------------------------------------------------------
// Newmark velocity and acceleration update
//   a^{n+1} = coeff * (u^{n+1} - u_tilde)
//   v^{n+1} = v^n + dt*((1-gamma)*a^n + gamma*a^{n+1})
// -----------------------------------------------------------------------
void EulerBeamInextensibleIndex1::update_newmark_state(real_t dt,
                                                               real_t beta,
                                                               real_t gamma)
{
  const real_t     coeff   = 1.0 / (beta * dt * dt);
  const VectorXd   u_tilde = u_prev + dt * v_prev + dt * dt * (0.5 - beta) * a_prev;
  const VectorXd   a_new   = coeff * (u - u_tilde);
  const VectorXd   v_new   = v_prev + dt * ((1.0 - gamma) * a_prev + gamma * a_new);
  a_prev = a_new;
  v_prev = v_new;
}

// -----------------------------------------------------------------------
// Zero velocity and acceleration at kinematically constrained DOFs
// -----------------------------------------------------------------------
void EulerBeamInextensibleIndex1::apply_dynamic_state_boundary_conditions()
{
  for (size_t bi = 0; bi < 2; ++bi) {
    const EulerBeamBCType bctype = boundary_conditions.type[bi];
    if (bctype != simple_bc && bctype != clamped_bc) continue;

    const size_t ni = (boundary_conditions.end[bi] == left) ? 0 : nodes - 1;
    const size_t ix = offset_x + 2 * ni;
    const size_t iy = offset_y + 2 * ni;
    const size_t iz = offset_z + 2 * ni;

    v_prev(ix) = v_prev(iy) = v_prev(iz) = 0.0;
    a_prev(ix) = a_prev(iy) = a_prev(iz) = 0.0;

    if (bctype == clamped_bc) {
      const size_t sx = offset_x + 2 * ni + 1;
      const size_t sy = offset_y + 2 * ni + 1;
      const size_t sz = offset_z + 2 * ni + 1;
      v_prev(sx) = v_prev(sy) = v_prev(sz) = 0.0;
      a_prev(sx) = a_prev(sy) = a_prev(sz) = 0.0;
    }
  }
}

// -----------------------------------------------------------------------
// Position drift correction
//
// Iterates Newton steps projecting u back onto C(u)=0 via
// the minimum-norm correction:
//
//   Solve  4*(B*B^T)*mu = c
//   u    -= 2*B^T * mu
//
// where c_a = int eta_a*(||r'||^2-1)ds and B_pos = 2*B.
// -----------------------------------------------------------------------
void EulerBeamInextensibleIndex1::project_position_onto_constraint()
{
  for (size_t iter = 0; iter < max_projection_iter; ++iter) {
    const VectorXd c      = assemble_constraint_residual_vector();
    const real_t   c_norm = c.norm();

    ELFF_LOG("  position projection iter " << iter
             << ": ||C|| = " << c_norm);

    if (c_norm < tol_position_drift) break;

    // Reassemble B at current (corrected) u
    SparseMatrix<real_t> B;
    VectorXd             g_dummy(ndof_l);
    assemble_B_and_g(B, g_dummy);

    // B_pos = 2*B,  B_pos*B_pos^T = 4*B*B^T
    const SparseMatrix<real_t> BBT_sparse = B * B.transpose();
    MatrixXd                   BBT        = 4.0 * BBT_sparse.toDense();
    BBT.diagonal().array() += 1e-14; // regularisation

    LLT<MatrixXd> llt;
    llt.compute(BBT);
    if (llt.info() != Success) {
      ELFF_WARNING("project_position_onto_constraint: B*B^T factorization "
                   "failed at iter "
                   << iter);
      break;
    }

    // delta_u = -B_pos^T * mu = -2*B^T * mu
    const VectorXd mu = llt.solve(c);
    u -= 2.0 * (B.transpose() * mu);

    // Re-enforce displacement BCs after correction
    for (size_t bi = 0; bi < 2; ++bi) {
      const size_t          ni     = (boundary_conditions.end[bi] == left) ? 0 : nodes - 1;
      const EulerBeamBCType bctype = boundary_conditions.type[bi];
      const EulerBeamBCVals bcvals = boundary_conditions.vals[bi];

      if (bctype == simple_bc || bctype == clamped_bc) {
        u(offset_x + 2 * ni) = bcvals.position[0];
        u(offset_y + 2 * ni) = bcvals.position[1];
        u(offset_z + 2 * ni) = bcvals.position[2];
      }
      if (bctype == clamped_bc) {
        u(offset_x + 2 * ni + 1) = bcvals.slope[0];
        u(offset_y + 2 * ni + 1) = bcvals.slope[1];
        u(offset_z + 2 * ni + 1) = bcvals.slope[2];
      }
    }
  }
}

// -----------------------------------------------------------------------
// Velocity drift correction
//
// Projects v_prev onto the velocity constraint manifold r' . r_dot' = 0:
//
//   Solve  B*B^T * mu = B*v_prev
//   v_prev -= B^T * mu
// -----------------------------------------------------------------------
void EulerBeamInextensibleIndex1::project_velocity_onto_constraint(
  const SparseMatrix<real_t>& B)
{
  const VectorXd dot_c      = B * v_prev;
  const real_t   drift_norm = dot_c.norm();

  if (drift_norm < tol_velocity_drift) return;

  ELFF_LOG("  velocity projection: ||B*v|| = " << drift_norm);

  const SparseMatrix<real_t> BBT_sparse = B * B.transpose();
  MatrixXd                   BBT        = BBT_sparse.toDense();
  BBT.diagonal().array() += 1e-14;

  LLT<MatrixXd> llt;
  llt.compute(BBT);
  if (llt.info() != Success) {
    ELFF_WARNING("project_velocity_onto_constraint: B*B^T factorization failed");
    return;
  }

  const VectorXd mu = llt.solve(dot_c);
  v_prev -= B.transpose() * mu;

  apply_dynamic_state_boundary_conditions();
}

// -----------------------------------------------------------------------
// Mesh update
// -----------------------------------------------------------------------
void EulerBeamInextensibleIndex1::update_mesh()
{
  std::vector<std::array<real_t, 3>>& centerline = mesh.get_centerline();
  std::vector<std::array<real_t, 3>>& slope = mesh.get_slope();
  std::vector<std::array<real_t, 3>>& vel = mesh.get_centerline_velocity();
  for (size_t i = 0; i < nodes; ++i) {
    centerline[i][0] = u(offset_x + 2 * i);
    centerline[i][1] = u(offset_y + 2 * i);
    centerline[i][2] = u(offset_z + 2 * i);
    slope[i][0] = u(offset_x + 2 * i + 1);
    slope[i][1] = u(offset_y + 2 * i + 1);
    slope[i][2] = u(offset_z + 2 * i + 1);
    vel[i][0] = v_prev(offset_x + 2 * i);
    vel[i][1] = v_prev(offset_y + 2 * i);
    vel[i][2] = v_prev(offset_z + 2 * i);
  }
}

} // namespace Models
} // namespace ELFF
