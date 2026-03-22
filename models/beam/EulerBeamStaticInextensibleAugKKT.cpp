#include "models/beam/EulerBeamStaticInextensibleAugKKT.hpp"

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
  , ndof(ndof_x + ndof_y + ndof_z + ndof_l)
  , offset_x(0)
  , offset_y(ndof_x)
  , offset_z(ndof_x + ndof_y)
  , offset_l(ndof_x + ndof_y + ndof_z)
  , r_penalty(r_penalty)
  , max_iter(1000)
  , tol(1e-6)
  , residual(VectorXd::Zero(ndof))
  , jacobian(MatrixXd::Zero(ndof, ndof))
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
  , ndof(ndof_x + ndof_y + ndof_z + ndof_l)
  , offset_x(0)
  , offset_y(ndof_x)
  , offset_z(ndof_x + ndof_y)
  , offset_l(ndof_x + ndof_y + ndof_z)
  , r_penalty(r_penalty)
  , max_iter(1000)
  , tol(1e-5)
  , residual(VectorXd::Zero(ndof))
  , jacobian(MatrixXd::Zero(ndof, ndof))
  , u(VectorXd::Zero(ndof))
{
  apply_initial_condition(this->mesh);
}

EulerBeamStaticInextensibleAugKKT::~EulerBeamStaticInextensibleAugKKT() =
  default;

void
EulerBeamStaticInextensibleAugKKT::solve()
{
  solve({ 0., 0., 0. });
}

void
EulerBeamStaticInextensibleAugKKT::solve(std::array<real_t, 3> load)
{
  for (size_t it = 0; it <= max_iter; it++) {
    assemble_system(load);
    apply_boundary_conditions();

    const real_t res_norm = residual.norm();

    if (res_norm < tol) {
      std::cout << "Converged in " << it << " iters.\n";
      break;
    } else if (it == max_iter) {
      ELFF_ABORT(
        "EulerBeamStaticInextensibleAugKKT::solve() did not converge.\n");
    }

    const VectorXd delta_u = jacobian.colPivHouseholderQr().solve(-residual);
    u += delta_u;
  }

  update_mesh();
}

void
EulerBeamStaticInextensibleAugKKT::apply_initial_condition()
{
  const size_t nodes = mesh.get_nodes();
  const real_t ds = mesh.get_ds();

  const size_t ndof_x = 2 * nodes;
  const size_t ndof_y = 2 * nodes;
  const size_t ndof_z = 2 * nodes;
  const size_t ndof_l = nodes;

  const size_t offset_x = 0;
  const size_t offset_y = ndof_x;
  const size_t offset_z = ndof_x + ndof_y;
  const size_t offset_l = ndof_x + ndof_y + ndof_z;

  for (size_t i = 0; i < nodes; i++) {
    u(offset_x + 2 * i + 0) = ds * i;
    u(offset_x + 2 * i + 1) = 1.;
    u(offset_y + 2 * i + 0) = 0.;
    u(offset_y + 2 * i + 1) = 0.;
    u(offset_z + 2 * i + 0) = 0.;
    u(offset_z + 2 * i + 1) = 0.;
    u(offset_l + i) = 0.;
  }
}

void
EulerBeamStaticInextensibleAugKKT::apply_initial_condition(EulerBeamMesh& bmesh)
{
  if (bmesh.get_nodes() != mesh.get_nodes()) {
    ELFF_ABORT("");
  }

  const auto centerline = bmesh.get_centerline();
  const auto slope = bmesh.get_slope();

  for (size_t i = 0; i < nodes; i++) {
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
EulerBeamStaticInextensibleAugKKT::assemble_residual(std::array<real_t, 3> load)
{
  residual = assemble_residual_template<real_t>(u, load);
}

void
EulerBeamStaticInextensibleAugKKT::assemble_system(std::array<real_t, 3> load)
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
EulerBeamStaticInextensibleAugKKT::apply_boundary_conditions()
{
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
EulerBeamStaticInextensibleAugKKT::update_mesh()
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
