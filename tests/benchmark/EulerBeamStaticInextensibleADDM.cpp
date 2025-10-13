

#include <gtest/gtest.h>
#include "models/beam/EulerBeamStaticInextensibleADDM.hpp"

using beam::EulerBeamStaticInextensibleADDM;
using beam::EulerBeam;
using beam::EulerBeamBCs;
using beam::real_t;
using beam::left;
using beam::right;
using beam::simple_bc;
using beam::clamped_bc;
using beam::free_bc;

static Eigen::Matrix4d K4_elem(double h) {
  Eigen::Matrix4d K;
  K <<  12.0/(h*h*h),  6.0/(h*h),   -12.0/(h*h*h),  6.0/(h*h),
        6.0/(h*h),     4.0/h,        -6.0/(h*h),    2.0/h,
       -12.0/(h*h*h), -6.0/(h*h),    12.0/(h*h*h), -6.0/(h*h),
        6.0/(h*h),     2.0/h,        -6.0/(h*h),    4.0/h;
  return K;
}

static Eigen::Matrix4d K2_elem(double h) {
  Eigen::Matrix4d K;
  K <<  6.0/(5.0*h),  1.0/10.0,   -6.0/(5.0*h),  1.0/10.0,
        1.0/10.0,     2.0*h/15.0, -1.0/10.0,    -h/30.0,
       -6.0/(5.0*h), -1.0/10.0,    6.0/(5.0*h), -1.0/10.0,
        1.0/10.0,    -h/30.0,     -1.0/10.0,     2.0*h/15.0;
  return K;
}

static Eigen::Matrix4d A_elem(double h, double EI, double R) {
  return EI * K4_elem(h) + R * K2_elem(h);
}

// ---- Assemble expected global A for one field (dense) on a uniform 1D mesh ----
static Eigen::MatrixXd assemble_expected_A_dense(std::size_t Nnodes,
                                                 double L, double EI, double R)
{
  // Nnodes = number of beam nodes; ndof = 2 per node (value, slope)
  const std::size_t ndof = 2 * Nnodes;
  const std::size_t Ne = Nnodes - 1;
  const double h = L / Ne;

  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(ndof, ndof);

  for (std::size_t e = 0; e < Ne; ++e) {
    // local dof map: [u_i, u'_i, u_{i+1}, u'_{i+1}]
    const int i0 = static_cast<int>(e);
    const int i1 = static_cast<int>(e + 1);
    const int edofs[4] = { 2*i0+0, 2*i0+1, 2*i1+0, 2*i1+1 };

    Eigen::Matrix4d Ae = A_elem(h, EI, R);

    for (int a = 0; a < 4; ++a) {
      int ra = edofs[a];
      for (int b = 0; b < 4; ++b) {
        int cb = edofs[b];
        A(ra, cb) += Ae(a, b);
      }
    }
  }
  return A;
}

// ---- Your production code should expose a function that returns A (dense) BEFORE BCs ----
// For the test, we assume you have something like:
//   Eigen::MatrixXd assemble_A_dense_raw(std::size_t Nnodes, double L, double EI, double R);
// which builds K4, K2 by quadrature and returns EI*K4 + R*K2 without touching BCs.

// If your current function applies BCs inside, consider splitting:
//   build_K4K2_dense(...) -> returns K4,K2
//   form_A_dense(K4,K2,EI,R) -> returns A
//   apply_BC_pattern_in_place(A, bc)
// so the test can compare the raw A.

TEST(EulerBeamStaticInextensibleADDMTest, MatchesClosedFormForSmallMesh)
{
 EulerBeamBCs boundary_conditions = {
    .end = {left, right},
    .type = {clamped_bc, free_bc},
    .vals = {{
      .position = {0,0,0},
      .slope = {1,0,0}
    }, {
      .position = {1,0,0},
      .slope = {1,0,0}
    }}
  };
  
  const std::size_t Nnodes = 3;   // 3 nodes -> 2 elements -> tiny, easy check
  const double      L      = 2.0; // so each element has h = 1.0
  const double      EI     = 3.7; // arbitrary nonzero
  const double      R      = 5.0; // penalty

  // Expected closed-form A (global, one field)
  Eigen::MatrixXd A_expected = assemble_expected_A_dense(Nnodes, L, EI, R);

  // Compute A via your numerical assembly (quadrature), **no BCs applied**
  EulerBeamStaticInextensibleADDM static_beam(L, EI, 0.0, 1.0, Nnodes, boundary_conditions, R);
  static_beam.assemble_A();
  Eigen::MatrixXd A_numerical = static_beam.get_A();

  // Symmetry sanity (optional)
  ASSERT_NEAR((A_numerical - A_numerical.transpose()).norm(), 0.0, 1e-12);
  ASSERT_NEAR((A_expected  - A_expected .transpose()).norm(),  0.0, 1e-12);

  // Entry-wise comparison
  const double tol = 1e-6;
  ASSERT_EQ(A_expected.rows(), A_numerical.rows());
  ASSERT_EQ(A_expected.cols(), A_numerical.cols());
  for (int i = 0; i < A_expected.rows(); ++i) {
    for (int j = 0; j < A_expected.cols(); ++j) {
      EXPECT_NEAR(A_expected(i,j), A_numerical(i,j), tol)
          << "Mismatch at (" << i << "," << j << ")";
    }
  }
}

TEST(EulerBeamStaticInextensibleADDMTest, SolveUniformLoadAndPlot) {

  real_t length = 1., EI = 1., area = 1., r_pentalty = 1e3, load = -10;
  size_t nodes = 15;

  EulerBeamBCs boundary_conditions = {
    .end = {left, right},
    .type = {clamped_bc, free_bc},
    .vals = {{
      .position = {0,0,0},
      .slope = {1,0,0}
    }, {
      .position = {1,0,0},
      .slope = {0,0,0}
    }}
  };

  EulerBeamStaticInextensibleADDM static_beam(length, EI, load, area, nodes, boundary_conditions, r_pentalty);

  //GTEST_LOG_(INFO) << "CTEST_FULL_OUTPUT";
  // static_beam.apply_initial_condition();
  // static_beam.update_mesh();
  // static_beam.plot("Static Inextensible Euler Beam");

  static_beam.solve();
  static_beam.plot("Static Inextensible Euler Beam");
  
  // auto u = static_beam.get_solution();
  // GTEST_LOG_(INFO) << "u =\n" << u << "\n";

  // Optionally, add some real checks:
  // EXPECT_NEAR(F(0), expected_value, 1e-12);
  // EXPECT_NEAR(K(0,0), expected_K00, 1e-12);
};