#include "models/beam/EulerBeamStaticInextensibleMoMSparse.hpp"

namespace ELFF {
namespace Models {

EulerBeamStaticInextensibleMoMSparse::EulerBeamStaticInextensibleMoMSparse(
  real_t length,
  real_t EI,
  size_t nodes,
  EulerBeam::EulerBeamBCs bcs,
  real_t r_penalty)
  : EulerBeam(length, EI, nodes, bcs)
  , dimension()
  , elements(nodes - 1)
  , nodes(nodes)
  , ds(mesh.get_ds())
  , ndof_x(2 * nodes)
  , ndof_y(2 * nodes)
  , ndof_z(2 * nodes)
  , ndof_l(1 * nodes)
  , offset_x(0)
  , offset_y(ndof_x)
  , offset_z(ndof_x + ndof_y)
  , offset_l()
  , ndof(ndof_x + ndof_y + ndof_z)
  , r_penalty(r_penalty)
  , max_iter_inner(1000)
  , max_iter_outer(1000)
  , tol_inner(1e-5)
  , tol_outer(1e-5)
  , residual(VectorXd::Zero(ndof))
  , lambda(VectorXd::Zero(ndof_l))
  , jacobian(SparseMatrix<real_t>(ndof, ndof))
  , u(VectorXd::Zero(ndof))
{
  apply_initial_condition(mesh);
}

EulerBeamStaticInextensibleMoMSparse::EulerBeamStaticInextensibleMoMSparse(
  real_t length,
  real_t EI,
  real_t mu,
  size_t nodes,
  EulerBeam::EulerBeamBCs bcs,
  real_t r_penalty)
  : EulerBeam(length, EI, mu, nodes, bcs)
  , dimension()
  , elements(nodes - 1)
  , nodes(nodes)
  , ds(mesh.get_ds())
  , ndof_x(2 * nodes)
  , ndof_y(2 * nodes)
  , ndof_z(2 * nodes)
  , ndof_l(1 * nodes)
  , offset_x(0)
  , offset_y(ndof_x)
  , offset_z(ndof_x + ndof_y)
  , offset_l()
  , ndof(ndof_x + ndof_y + ndof_z)
  , r_penalty(r_penalty)
  , max_iter_inner(1000)
  , max_iter_outer(1000)
  , tol_inner(1e-6)
  , tol_outer(1e-6)
  , residual(VectorXd::Zero(ndof))
  , lambda(VectorXd::Zero(ndof_l))
  , jacobian(SparseMatrix<real_t>(ndof, ndof))
  , u(VectorXd::Zero(ndof))
{
  apply_initial_condition(mesh);
}

EulerBeamStaticInextensibleMoMSparse::~EulerBeamStaticInextensibleMoMSparse() =
  default;

void
EulerBeamStaticInextensibleMoMSparse::solve()
{
  solve({ 0., 0., 0. });
}

void
EulerBeamStaticInextensibleMoMSparse::solve(std::array<real_t, 3> load)
{
  real_t S_norm = 0;

  ConjugateGradient<SparseMatrix<real_t>,
                    Lower | Upper,
                    IncompleteCholesky<real_t>>
    solver;

  for (size_t iter_outer = 0; iter_outer < max_iter_outer; iter_outer++) {
    assemble_system(load);
    apply_boundary_conditions();

    const real_t res_norm = residual.norm();

    if (res_norm < tol_outer) {
      break;
    } else if (iter_outer == max_iter_outer - 1) {
      ELFF_ABORT("EulerBeamStaticInexntensibleMoM::solve() did not converge.\n");
    }

    solver.setTolerance(tol_inner);
    solver.compute(jacobian);

    if (solver.info() != Success) {
      ELFF_ABORT("EulerBeamStaticInextensibleMoMSparse::solve(): "
                 "Preconditioner failed.\n");
    }

    const VectorXd delta_u = solver.solve(-residual);
    u += delta_u;

    S_norm = update_lambda();
  }

  (void) S_norm;
  update_mesh();
}

void
EulerBeamStaticInextensibleMoMSparse::solve(
  std::vector<std::array<real_t, 3>> load)
{
  ELFF_ASSERT(load.size() == nodes,
              "Size of load vector must equal number of nodes.");

  real_t S_norm = 0;

  ConjugateGradient<SparseMatrix<real_t>,
                    Lower | Upper,
                    IncompleteCholesky<real_t>>
    solver;

  for (size_t iter_outer = 0; iter_outer < max_iter_outer; iter_outer++) {
    assemble_system(load);
    apply_boundary_conditions();

    const real_t res_norm = residual.norm();

    if (res_norm < tol_outer) {
      break;
    } else if (iter_outer == max_iter_outer - 1) {
      ELFF_ABORT("EulerBeamStaticInexntensibleMoM::solve() did not converge.\n");
    }

    solver.setTolerance(tol_inner);
    solver.compute(jacobian);

    if (solver.info() != Success) {
      ELFF_ABORT("EulerBeamStaticInextensibleMoMSparse::solve(): "
                 "Preconditioner failed.\n");
    }

    const VectorXd delta_u = solver.solve(-residual);
    u += delta_u;

    S_norm = update_lambda();
  }

  (void) S_norm;
  update_mesh();
}

void
EulerBeamStaticInextensibleMoMSparse::apply_initial_condition(
  EulerBeamMesh& bmesh)
{
  ELFF_ASSERT(
    nodes == bmesh.get_nodes(),
    "Provided mesh must have same node count as the previous mesh.\n");

  const auto centerline = bmesh.get_centerline();
  const auto slopes = bmesh.get_slope();

  for (size_t ni = 0; ni < nodes; ni++) {
    u(offset_x + 2 * ni + 0) = centerline[ni][0];
    u(offset_x + 2 * ni + 1) = slopes[ni][0];
    u(offset_y + 2 * ni + 0) = centerline[ni][1];
    u(offset_y + 2 * ni + 1) = slopes[ni][1];
    u(offset_z + 2 * ni + 0) = centerline[ni][2];
    u(offset_z + 2 * ni + 1) = slopes[ni][2];
  }

  update_mesh();
}

void
EulerBeamStaticInextensibleMoMSparse::apply_initial_condition()
{
  for (size_t i = 0; i < nodes; i++) {
    u(offset_x + 2 * i + 0) = ds * i;
    u(offset_x + 2 * i + 1) = 1.;
    u(offset_y + 2 * i + 0) = 0.;
    u(offset_y + 2 * i + 1) = 0.;
    u(offset_z + 2 * i + 0) = 0.;
    u(offset_z + 2 * i + 1) = 0.;
  }
}

void
EulerBeamStaticInextensibleMoMSparse::assemble_residual(
  std::array<real_t, 3> load)
{
  residual = assemble_residual_template<real_t>(u, load);
}

void
EulerBeamStaticInextensibleMoMSparse::assemble_residual(
  std::vector<std::array<real_t, 3>> load)
{
  residual = assemble_residual_template<real_t>(u, load);
}

std::array<size_t, 12>
EulerBeamStaticInextensibleMoMSparse::get_element_dof_indices(size_t e) const
{
  const size_t n0 = e;
  const size_t n1 = e + 1;

  return { offset_x + 2 * n0 + 0,
           offset_x + 2 * n0 + 1,
           offset_x + 2 * n1 + 0,
           offset_x + 2 * n1 + 1,
           offset_y + 2 * n0 + 0,
           offset_y + 2 * n0 + 1,
           offset_y + 2 * n1 + 0,
           offset_y + 2 * n1 + 1,
           offset_z + 2 * n0 + 0,
           offset_z + 2 * n0 + 1,
           offset_z + 2 * n1 + 0,
           offset_z + 2 * n1 + 1 };
}

Matrix<real_t, 12, 1>
EulerBeamStaticInextensibleMoMSparse::get_element_state(
  const std::array<size_t, 12>& idx) const
{
  Matrix<real_t, 12, 1> u_elem;
  for (int i = 0; i < 12; ++i) {
    u_elem(i) = u(idx[i]);
  }
  return u_elem;
}

std::array<real_t, 2>
EulerBeamStaticInextensibleMoMSparse::get_element_lambda(size_t e) const
{
  return { lambda(e), lambda(e + 1) };
}

real_t
EulerBeamStaticInextensibleMoMSparse::update_lambda(real_t omega)
{
  (void) omega;

  const real_t xi_q[] = { 0.1127016654, 0.5, 0.8872983346 };
  const real_t w_q[] = { 0.2777777778, 0.4444444444, 0.2777777778 };

  Matrix<real_t, Dynamic, 1> lambda_n =
    Matrix<real_t, Dynamic, 1>::Zero(ndof_l);

  for (size_t e = 0; e < elements; ++e) {
    const std::vector<size_t> elem_nodes = { e, e + 1 };
    const std::vector<size_t> idx_x = { offset_x + 2 * elem_nodes[0],
                                        offset_x + 2 * elem_nodes[0] + 1,
                                        offset_x + 2 * elem_nodes[1],
                                        offset_x + 2 * elem_nodes[1] + 1 };
    const std::vector<size_t> idx_y = { offset_y + 2 * elem_nodes[0],
                                        offset_y + 2 * elem_nodes[0] + 1,
                                        offset_y + 2 * elem_nodes[1],
                                        offset_y + 2 * elem_nodes[1] + 1 };
    const std::vector<size_t> idx_z = { offset_z + 2 * elem_nodes[0],
                                        offset_z + 2 * elem_nodes[0] + 1,
                                        offset_z + 2 * elem_nodes[1],
                                        offset_z + 2 * elem_nodes[1] + 1 };
    const std::vector<size_t> idx_l = { elem_nodes[0], elem_nodes[1] };

    const std::array<real_t, 4> ux = {
      u[idx_x[0]], u[idx_x[1]], u[idx_x[2]], u[idx_x[3]]
    };
    const std::array<real_t, 4> uy = {
      u[idx_y[0]], u[idx_y[1]], u[idx_y[2]], u[idx_y[3]]
    };
    const std::array<real_t, 4> uz = {
      u[idx_z[0]], u[idx_z[1]], u[idx_z[2]], u[idx_z[3]]
    };
    const std::array<real_t, 2> ul = { lambda[idx_l[0]], lambda[idx_l[1]] };

    std::vector<real_t> R_loc_l(2, 0);

    for (size_t qi = 0; qi < 3; ++qi) {
      const real_t xi = xi_q[qi];
      const real_t w = w_q[qi];

      const auto H = CubicHermite<real_t>::values(xi, ds);
      const auto dH = CubicHermite<real_t>::derivs(xi, ds);
      const auto ddH = CubicHermite<real_t>::second_derivs(xi, ds);
      const auto M = LinearShape<real_t>::values(xi);

      real_t x = 0, xp = 0, xpp = 0;
      real_t y = 0, yp = 0, ypp = 0;
      real_t z = 0, zp = 0, zpp = 0;

      for (size_t i = 0; i < 4; i++) {
        x += H[i] * ux[i];
        xp += dH[i] * ux[i];
        xpp += ddH[i] * ux[i];
        y += H[i] * uy[i];
        yp += dH[i] * uy[i];
        ypp += ddH[i] * uy[i];
        z += H[i] * uz[i];
        zp += dH[i] * uz[i];
        zpp += ddH[i] * uz[i];
      }

      real_t l = 0;
      for (size_t i = 0; i < 2; i++) {
        l += M[i] * ul[i];
      }

      const real_t S = xp * xp + yp * yp + zp * zp - 1.0;

      for (size_t a = 0; a < 2; ++a) {
        R_loc_l[a] += S * M[a] * w * ds;
      }
    }

    for (size_t i = 0; i < 2; ++i) {
      lambda_n[idx_l[i]] += R_loc_l[i];
    }
  }

  lambda = lambda_n;
  return lambda_n.norm();
}

void
EulerBeamStaticInextensibleMoMSparse::update_mesh()
{
  const size_t nodes = this->mesh.get_nodes();

  std::vector<std::array<real_t, 3>>& centerline = mesh.get_centerline();
  std::vector<std::array<real_t, 3>>& slope = mesh.get_slope();
  std::vector<real_t>& s = mesh.get_curvilinear_axis();
  (void)s;

  for (size_t i = 0; i < nodes; ++i) {
    centerline[i][0] = u(offset_x + 2 * i);
    centerline[i][1] = u(offset_y + 2 * i);
    centerline[i][2] = u(offset_z + 2 * i);
    slope[i][0] = u(offset_x + 2 * i + 1);
    slope[i][1] = u(offset_y + 2 * i + 1);
    slope[i][2] = u(offset_z + 2 * i + 1);
  }
}

void
EulerBeamStaticInextensibleMoMSparse::assemble_system(
  std::array<real_t, 3> load)
{
  using ADDeriv = Matrix<real_t, 12, 1>;
  using AD = AutoDiffScalar<ADDeriv>;
  using ADVec = Matrix<AD, 12, 1>;
  using Tpl = Triplet<real_t>;

  residual = VectorXd::Zero(ndof);
  std::vector<Tpl> triplets;
  triplets.reserve(elements * 12 * 12);

  for (size_t e = 0; e < elements; ++e) {
    const auto idx = get_element_dof_indices(e);
    const auto lambda_elem = get_element_lambda(e);
    const Matrix<real_t, 12, 1> u_elem = get_element_state(idx);

    ADVec u_ad;
    for (int a = 0; a < 12; ++a) {
      ADDeriv seed = ADDeriv::Zero();
      seed(a) = 1.0;
      u_ad(a) = AD(u_elem(a), seed);
    }

    const ADVec R_loc_ad =
      assemble_element_residual_template<AD>(u_ad, lambda_elem, load);

    for (int a = 0; a < 12; ++a) {
      residual(idx[a]) += R_loc_ad(a).value();

      const ADDeriv& dRa = R_loc_ad(a).derivatives();
      for (int b = 0; b < 12; ++b) {
        const real_t dj = dRa(b);
        if (dj != 0.0) {
          triplets.emplace_back(idx[a], idx[b], dj);
        }
      }
    }
  }

  jacobian.resize(ndof, ndof);
  jacobian.setFromTriplets(triplets.begin(), triplets.end());
  jacobian.makeCompressed();
}

void
EulerBeamStaticInextensibleMoMSparse::assemble_system(
  std::vector<std::array<real_t, 3>> load)
{
  using ADDeriv = Matrix<real_t, 12, 1>;
  using AD = AutoDiffScalar<ADDeriv>;
  using ADVec = Matrix<AD, 12, 1>;
  using Tpl = Triplet<real_t>;

  residual = VectorXd::Zero(ndof);
  std::vector<Tpl> triplets;
  triplets.reserve(elements * 12 * 12);

  for (size_t e = 0; e < elements; ++e) {
    const auto idx = get_element_dof_indices(e);
    const auto lambda_elem = get_element_lambda(e);
    const Matrix<real_t, 12, 1> u_elem = get_element_state(idx);
    const std::array<std::array<real_t, 3>, 2> load_elem = { load[e], load[e + 1] };

    ADVec u_ad;
    for (int a = 0; a < 12; ++a) {
      ADDeriv seed = ADDeriv::Zero();
      seed(a) = 1.0;
      u_ad(a) = AD(u_elem(a), seed);
    }

    const ADVec R_loc_ad =
      assemble_element_residual_template<AD>(u_ad, lambda_elem, load_elem);

    for (int a = 0; a < 12; ++a) {
      residual(idx[a]) += R_loc_ad(a).value();

      const ADDeriv& dRa = R_loc_ad(a).derivatives();
      for (int b = 0; b < 12; ++b) {
        const real_t dj = dRa(b);
        if (dj != 0.0) {
          triplets.emplace_back(idx[a], idx[b], dj);
        }
      }
    }
  }

  jacobian.resize(ndof, ndof);
  jacobian.setFromTriplets(triplets.begin(), triplets.end());
  jacobian.makeCompressed();
}

void
EulerBeamStaticInextensibleMoMSparse::apply_boundary_conditions()
{
  for (size_t bi = 0; bi < 2; ++bi) {
    const EulerBeamBCEnd bcend = boundary_conditions.end[bi];
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

    const EulerBeamBCType bctype = boundary_conditions.type[bi];
    const EulerBeamBCVals bcvals = boundary_conditions.vals[bi];
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
        for (size_t i = 0; i < vals.size(); i++) {
          residual[idx[i]] -= vals[i];
        }
        break;
      default:
        for (auto i : idx) {
          for (SparseMatrix<real_t>::InnerIterator it(jacobian, i); it; ++it) {
            it.valueRef() = 0.0;
          }

          for (int col = 0; col < jacobian.outerSize(); ++col) {
            for (SparseMatrix<real_t>::InnerIterator it(jacobian, col); it;
                 ++it) {
              if (it.row() == i) {
                it.valueRef() = 0.0;
              }
            }
          }

          jacobian.coeffRef(i, i) = 1.0;
        }

        for (size_t i = 0; i < vals.size(); i++) {
          residual[idx[i]] = u[idx[i]] - vals[i];
        }
        break;
    }
  }
}

} // namespace Models
} // namespace ELFF
