#include "elff/models/beam/EulerBeamStaticInextensibleAugKKT.hpp"

namespace ELFF {
namespace Models {

EulerBeamStaticInextensibleAugKKT::EulerBeamStaticInextensibleAugKKT(
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
  , offset_x(0)
  , offset_y(ndof_x)
  , offset_z(ndof_x + ndof_y)
  , offset_l(ndof_x + ndof_y + ndof_z)
  , ndof(ndof_x + ndof_y + ndof_z + ndof_l)
  , r_penalty(r_penalty)
  , max_iter(20)
  , tol(1e-5)
  , residual(VectorXd::Zero(ndof))
  , jacobian(SparseMatrix<real_t>(ndof, ndof))
  , u(VectorXd::Zero(ndof))
{
  apply_initial_condition(this->mesh);
}

EulerBeamStaticInextensibleAugKKT::EulerBeamStaticInextensibleAugKKT(
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
  , offset_x(0)
  , offset_y(ndof_x)
  , offset_z(ndof_x + ndof_y)
  , offset_l(ndof_x + ndof_y + ndof_z)
  , ndof(ndof_x + ndof_y + ndof_z + ndof_l)
  , r_penalty(r_penalty)
  , max_iter(100)
  , tol(1e-5)
  , residual(VectorXd::Zero(ndof))
  , jacobian(SparseMatrix<real_t>(ndof, ndof))
  , u(VectorXd::Zero(ndof))
{
  apply_initial_condition(this->mesh);
}

EulerBeamStaticInextensibleAugKKT::
  ~EulerBeamStaticInextensibleAugKKT() = default;

void
EulerBeamStaticInextensibleAugKKT::solve()
{
  solve({ 0., 0., 0. });
}

void
EulerBeamStaticInextensibleAugKKT::solve(std::array<real_t, 3> load)
{
  SparseLU<SparseMatrix<real_t>, COLAMDOrdering<int>> solver;

  for (size_t it = 0; it <= max_iter; ++it) {
    assemble_system(load);
    apply_boundary_conditions();

    const real_t res_norm = residual.norm();

    if (res_norm < tol) {
      ELFF_LOG("Converged in " << it << " iters.");
      break;
    } else if (it == max_iter) {
      ELFF_ABORT(
        "EulerBeamStaticInextensibleAugKKT::solve() did not converge.\n");
    }

    solver.analyzePattern(jacobian);
    solver.factorize(jacobian);

    if (solver.info() != Success) {
      ELFF_ABORT("EulerBeamStaticInextensibleAugKKT::solve(): "
                 "SparseLU factorization failed.\n");
    }

    const VectorXd delta_u = solver.solve(-residual);

    if (solver.info() != Success) {
      ELFF_ABORT("EulerBeamStaticInextensibleAugKKT::solve(): "
                 "SparseLU solve failed.\n");
    }

    u += delta_u;
  }

  update_mesh();
}

void
EulerBeamStaticInextensibleAugKKT::apply_initial_condition()
{
  for (size_t i = 0; i < nodes; ++i) {
    u(offset_x + 2 * i + 0) = ds * i;
    u(offset_x + 2 * i + 1) = 1.;
    u(offset_y + 2 * i + 0) = 0.;
    u(offset_y + 2 * i + 1) = 0.;
    u(offset_z + 2 * i + 0) = 0.;
    u(offset_z + 2 * i + 1) = 0.;
    u(offset_l + i) = 0.;
  }

  update_mesh();
}

void
EulerBeamStaticInextensibleAugKKT::apply_initial_condition(
  EulerBeamMesh& bmesh)
{
  if (bmesh.get_nodes() != mesh.get_nodes()) {
    ELFF_ABORT("");
  }

  const auto centerline = bmesh.get_centerline();
  const auto slope = bmesh.get_slope();

  for (size_t i = 0; i < nodes; ++i) {
    u(offset_x + 2 * i + 0) = centerline[i][0];
    u(offset_x + 2 * i + 1) = slope[i][0];
    u(offset_y + 2 * i + 0) = centerline[i][1];
    u(offset_y + 2 * i + 1) = slope[i][1];
    u(offset_z + 2 * i + 0) = centerline[i][2];
    u(offset_z + 2 * i + 1) = slope[i][2];
    u(offset_l + i) = 0.;
  }
}

void
EulerBeamStaticInextensibleAugKKT::assemble_residual(
  std::array<real_t, 3> load)
{
  residual = VectorXd::Zero(ndof);

  for (size_t e = 0; e < elements; ++e) {
    const auto idx = get_element_dof_indices(e);
    const auto u_elem = get_element_state(idx);
    const auto r_elem = assemble_element_residual_template<real_t>(u_elem, load);

    for (int a = 0; a < 14; ++a) {
      residual(idx[a]) += r_elem(a);
    }
  }
}

std::array<size_t, 14>
EulerBeamStaticInextensibleAugKKT::get_element_dof_indices(size_t e) const
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
           offset_z + 2 * n1 + 1,
           offset_l + n0,
           offset_l + n1 };
}

Matrix<real_t, 14, 1>
EulerBeamStaticInextensibleAugKKT::get_element_state(
  const std::array<size_t, 14>& idx) const
{
  Matrix<real_t, 14, 1> u_elem;
  for (int i = 0; i < 14; ++i) {
    u_elem(i) = u(idx[i]);
  }
  return u_elem;
}

void
EulerBeamStaticInextensibleAugKKT::assemble_system(
  std::array<real_t, 3> load)
{
  using ADDeriv = Matrix<real_t, 14, 1>;
  using AD = AutoDiffScalar<ADDeriv>;
  using ADVec = Matrix<AD, 14, 1>;
  using Tpl = Triplet<real_t>;

  residual = VectorXd::Zero(ndof);
  std::vector<Tpl> triplets;
  triplets.reserve(elements * 14 * 14);

  for (size_t e = 0; e < elements; ++e) {
    const auto idx = get_element_dof_indices(e);
    const Matrix<real_t, 14, 1> u_elem = get_element_state(idx);

    ADVec u_ad;
    for (int a = 0; a < 14; ++a) {
      ADDeriv seed = ADDeriv::Zero();
      seed(a) = 1.0;
      u_ad(a) = AD(u_elem(a), seed);
    }

    const ADVec r_elem_ad = assemble_element_residual_template<AD>(u_ad, load);

    for (int a = 0; a < 14; ++a) {
      residual(idx[a]) += r_elem_ad(a).value();

      const ADDeriv& dRa = r_elem_ad(a).derivatives();
      for (int b = 0; b < 14; ++b) {
        const real_t value = dRa(b);
        if (value != 0.0) {
          triplets.emplace_back(idx[a], idx[b], value);
        }
      }
    }
  }

  jacobian.resize(ndof, ndof);
  jacobian.setFromTriplets(triplets.begin(), triplets.end());
  jacobian.makeCompressed();
}

void
EulerBeamStaticInextensibleAugKKT::apply_boundary_conditions()
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
        idx = { offset_x + 2 * ni + 0,
                offset_y + 2 * ni + 0,
                offset_z + 2 * ni + 0 };
        vals = { bcvals.force[0], bcvals.force[1], bcvals.force[2] };
        break;
      case point_torque_bc:
        idx = { offset_x + 2 * ni + 1,
                offset_y + 2 * ni + 1,
                offset_z + 2 * ni + 1 };
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
        for (auto i : idx) {
          for (SparseMatrix<real_t>::InnerIterator it(jacobian, i); it; ++it) {
            it.valueRef() = 0.0;
          }

          for (int col = 0; col < jacobian.outerSize(); ++col) {
            for (SparseMatrix<real_t>::InnerIterator it(jacobian, col); it;
                 ++it) {
              if (it.row() == static_cast<int>(i)) {
                it.valueRef() = 0.0;
              }
            }
          }

          jacobian.coeffRef(i, i) = 1.0;
        }

        for (size_t i = 0; i < vals.size(); ++i) {
          residual(idx[i]) = u[idx[i]] - vals[i];
        }
        break;
    }
  }

  for (size_t bi = 0; bi < 2; ++bi) {
    const EulerBeamBCType bctype = boundary_conditions.type[bi];
    if (bctype != free_bc && bctype != simple_bc && bctype != clamped_bc) {
      continue;
    }

    size_t li = offset_l;
    switch (boundary_conditions.end[bi]) {
      case left:
        li += 0;
        break;
      case right:
        li += nodes - 1;
        break;
    }

    for (SparseMatrix<real_t>::InnerIterator it(jacobian, li); it; ++it) {
      it.valueRef() = 0.0;
    }

    for (int col = 0; col < jacobian.outerSize(); ++col) {
      for (SparseMatrix<real_t>::InnerIterator it(jacobian, col); it; ++it) {
        if (it.row() == static_cast<int>(li)) {
          it.valueRef() = 0.0;
        }
      }
    }

    jacobian.coeffRef(li, li) = 1.0;
    residual(li) = u(li);
  }
}

void
EulerBeamStaticInextensibleAugKKT::update_mesh()
{
  std::vector<std::array<real_t, 3>>& centerline = mesh.get_centerline();
  std::vector<std::array<real_t, 3>>& slope = mesh.get_slope();

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
