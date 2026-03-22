#pragma once

#include <Eigen/Dense>
#include <array>
#include <vector>
#include <utility>
#include <cmath>

namespace ELFF {
namespace Models {

using Vec3 = Eigen::Vector3d;
using Vec2 = Eigen::Vector2d;
using Mat2 = Eigen::Matrix2d;
using Mat32 = Eigen::Matrix<double, 3, 2>;
using Mat23 = Eigen::Matrix<double, 2, 3>;
using Mat3 = Eigen::Matrix3d;

// -----------------------------------------------------------------------------
// Domain / periodic helpers
// -----------------------------------------------------------------------------

struct PeriodicBox {
  Vec3 origin = Vec3::Zero();
  Vec3 length = Vec3::Ones();
  std::array<bool,3> periodic{true, true, true};

  Vec3 shortestDelta(const Vec3& a, const Vec3& b) const {
    Vec3 d = a - b;
    for (int k = 0; k < 3; ++k) {
      if (!periodic[k]) continue;
      const double L = length[k];
      if (d[k] >  0.5*L) d[k] -= L;
      if (d[k] < -0.5*L) d[k] += L;
    }
    return d;
  }

  Vec3 wrap(const Vec3& x) const {
    Vec3 y = x;
    for (int k = 0; k < 3; ++k) {
      if (!periodic[k]) continue;
      const double xmin = origin[k];
      const double xmax = origin[k] + length[k];
      while (y[k] >= xmax) y[k] -= length[k];
      while (y[k] <  xmin) y[k] += length[k];
    }
    return y;
  }
};

// -----------------------------------------------------------------------------
// Topology
// -----------------------------------------------------------------------------

struct NodeTopology {
  std::vector<int> neighbors;   // 5 or 6
  std::vector<int> edges;
  std::vector<int> triangles;
};

struct Edge {
  std::array<int,2> nodes{};
  std::array<int,2> triangles{};
  double restLength = 0.0;
};

struct TriangleTopology {
  std::array<int,3> nodes{};
  std::array<int,3> edges{};
};

// -----------------------------------------------------------------------------
// Reference data stored once after initialization
// -----------------------------------------------------------------------------

struct TriangleReference {
  Mat2 refShape = Mat2::Zero();          // columns = node1,node2 in reference plane
  Eigen::Matrix<double,3,2> dNdXi = Eigen::Matrix<double,3,2>::Zero();
  double restArea = 0.0;
};

// -----------------------------------------------------------------------------
// Current geometry/state
// -----------------------------------------------------------------------------

struct TriangleGeometry {
  double area = 0.0;
  Vec3 normal = Vec3::Zero();
  Vec3 centroid = Vec3::Zero();
  Vec2 stretch = Vec2::Zero();
  Vec2 tension = Vec2::Zero();
};

struct MeshState {
  Eigen::Matrix3Xd x;        // node positions
  Eigen::Matrix3Xd v;        // node velocities
  Eigen::Matrix3Xd f;        // nodal forces
  Eigen::Matrix3Xd normal;   // nodal normals

  Eigen::VectorXd meanCurv;
  Eigen::VectorXd refCurv;
  Eigen::VectorXd gaussCurv;

  Vec3 centroid = Vec3::Zero();
  double volume = 0.0;
  double initialVolume = 0.0;

  std::vector<TriangleGeometry> triGeom;

  void resize(int nNodes, int nTriangles) {
    x.setZero(3, nNodes);
    v.setZero(3, nNodes);
    f.setZero(3, nNodes);
    normal.setZero(3, nNodes);
    meanCurv.setZero(nNodes);
    refCurv.setZero(nNodes);
    gaussCurv.setZero(nNodes);
    triGeom.resize(nTriangles);
  }

  void zeroForces() { f.setZero(); }
};

struct CapsuleMesh3D {
  PeriodicBox box;
  std::vector<NodeTopology> nodeTopo;
  std::vector<Edge> edges;
  std::vector<TriangleTopology> triangles;
  std::vector<TriangleReference> triRef;
  MeshState state;

  int numNodes() const { return static_cast<int>(nodeTopo.size()); }
  int numEdges() const { return static_cast<int>(edges.size()); }
  int numTriangles() const { return static_cast<int>(triangles.size()); }
};

// -----------------------------------------------------------------------------
// Constitutive law concept
// -----------------------------------------------------------------------------

struct NeoHookeanLaw {
  double Es = 1.0;

  Vec2 dWdLambda(const Vec2& lambda) const {
    const double l1 = lambda[0], l2 = lambda[1];
    return {
      Es/(3.0*l1)*(l1*l1 - 1.0/std::pow(l1*l2, 2)),
      Es/(3.0*l2)*(l2*l2 - 1.0/std::pow(l1*l2, 2))
    };
  }
};

struct SkalakLaw {
  double Es = 1.0;
  double C = 1.0;

  Vec2 dWdLambda(const Vec2& lambda) const {
    const double l1 = lambda[0], l2 = lambda[1];
    const double J2m1 = std::pow(l1*l2, 2) - 1.0;
    return {
      Es*(l1*(l1*l1 - 1.0) + C*l1*l2*l2*J2m1),
      Es*(l2*(l2*l2 - 1.0) + C*l2*l1*l1*J2m1)
    };
  }
};

// -----------------------------------------------------------------------------
// Geometry ops
// -----------------------------------------------------------------------------

class GeometryOps3D {
public:
  static void updateTriangleGeometry(CapsuleMesh3D& mesh);
  static void updateNodeNormals(CapsuleMesh3D& mesh);
  static void updateCentroid(CapsuleMesh3D& mesh);
  static void updateVolume(CapsuleMesh3D& mesh);
  static double edgeLength(const CapsuleMesh3D& mesh, int eid);

  // current triangle -> local reference plane
  static std::pair<Mat32, Mat3> rotateTriangleToReferencePlane(
    const CapsuleMesh3D& mesh, int tid);
};

// -----------------------------------------------------------------------------
// Elasticity assembler
// -----------------------------------------------------------------------------

template<class Law>
class ElasticityAssembler3D {
public:
  explicit ElasticityAssembler3D(Law law) : law_(std::move(law)) {}

  void initializeReferenceConfiguration(CapsuleMesh3D& mesh) const {
    GeometryOps3D::updateTriangleGeometry(mesh);
    mesh.triRef.resize(mesh.numTriangles());

    for (int tid = 0; tid < mesh.numTriangles(); ++tid) {
      const auto [refNodes, RcurToRef] =
        GeometryOps3D::rotateTriangleToReferencePlane(mesh, tid);

      auto& ref = mesh.triRef[tid];
      ref.refShape.col(0) = refNodes.col(0).head<2>();
      ref.refShape.col(1) = refNodes.col(1).head<2>();
      ref.restArea = mesh.state.triGeom[tid].area;

      const double x1 = ref.refShape(0,0), y1 = ref.refShape(1,0);
      const double x2 = ref.refShape(0,1), y2 = ref.refShape(1,1);
      const double det = x1*y2 - y1*x2;

      ref.dNdXi.row(0) << (y1 - y2)/det, (x2 - x1)/det;
      ref.dNdXi.row(1) << y2/det,        -x2/det;
      ref.dNdXi.row(2) << -y1/det,        x1/det;
    }
  }

  void assemble(CapsuleMesh3D& mesh) const {
    GeometryOps3D::updateTriangleGeometry(mesh);
    GeometryOps3D::updateNodeNormals(mesh);
    mesh.state.zeroForces();

    for (int tid = 0; tid < mesh.numTriangles(); ++tid) {
      const auto& tri = mesh.triangles[tid];
      const auto& ref = mesh.triRef[tid];
      auto& geom = mesh.state.triGeom[tid];

      const auto [curNodesInRefPlane, RrefToCur] =
        GeometryOps3D::rotateTriangleToReferencePlane(mesh, tid);

      Mat2 U = Mat2::Zero();
      U.col(0) = curNodesInRefPlane.col(0).head<2>() - ref.refShape.col(0);
      U.col(1) = curNodesInRefPlane.col(1).head<2>() - ref.refShape.col(1);

      Mat2 F = Mat2::Identity();
      for (int a = 1; a < 3; ++a)
        F += U.col(a - 1) * ref.dNdXi.row(a);

      Mat2 C = F.transpose() * F;
      Eigen::SelfAdjointEigenSolver<Mat2> es(C);
      Vec2 eig = es.eigenvalues().cwiseMax(0.0);
      Vec2 lambda = eig.cwiseSqrt();

      Vec2 dW = law_.dWdLambda(lambda);
      geom.stretch = lambda;
      geom.tension << dW[0]/lambda[1], dW[1]/lambda[0];

      for (int a = 0; a < 3; ++a) {
        Vec2 fPlane = localNodalForce(a, F, C, lambda, ref.dNdXi, dW);
        Vec3 f3 = geom.area * RrefToCur.leftCols<2>() * fPlane;
        mesh.state.f.col(tri.nodes[a]) -= f3;
      }
    }
  }

private:
  Law law_;

  static Vec2 localNodalForce(
    int a,
    const Mat2& F,
    const Mat2& C,
    const Vec2& lambda,
    const Eigen::Matrix<double,3,2>& dNdXi,
    const Vec2& dW);
};

} // namespace model
} // namespace ELFF