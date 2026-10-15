#pragma once

#include <Eigen/Dense>
#include <array>
#include <vector>

namespace ELFF {
namespace Models {

using Vec3 = Eigen::Vector3d;
using Vec2 = Eigen::Vector2d;
using Mat2 = Eigen::Matrix2d;
using Mat32 = Eigen::Matrix<double, 3, 2>;
using Mat23 = Eigen::Matrix<double, 2, 3>;
using Mat3 = Eigen::Matrix3d;

struct CapsuleNodeTopology {
  std::vector<int> neighbors;
  std::vector<int> edges;
  std::vector<int> triangles;
};

struct CapsuleEdge {
  std::array<int,2> nodes{};
  std::array<int,2> triangles{};
  double restLength = 0.0;
};

struct CapsuleTriangleTopology {
  std::array<int,3> nodes{};
  std::array<int,3> edges{};
};

struct CapsuleTriangleReference {
  Mat2 refShape = Mat2::Zero();
  Eigen::Matrix<double,3,2> dNdXi = Eigen::Matrix<double,3,2>::Zero();
  double restArea = 0.0;
};

struct CapsuleTriangleGeometry {
  double area = 0.0;
  Vec3 normal = Vec3::Zero();
  Vec3 centroid = Vec3::Zero();
  Vec2 stretch = Vec2::Zero();
  Vec2 tension = Vec2::Zero();
};

struct CapsuleMeshState {
  Eigen::Matrix3Xd x;
  Eigen::Matrix3Xd v;
  Eigen::Matrix3Xd f;
  Eigen::Matrix3Xd normal;

  Eigen::VectorXd meanCurv;
  Eigen::VectorXd refCurv;
  Eigen::VectorXd gaussCurv;

  Vec3 centroid = Vec3::Zero();
  double volume = 0.0;
  double initialVolume = 0.0;

  std::vector<CapsuleTriangleGeometry> triGeom;

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

  void zeroForces() {
    f.setZero();
  }
};

struct CapsuleMesh {
  std::vector<CapsuleNodeTopology> nodeTopo;
  std::vector<CapsuleEdge> edges;
  std::vector<CapsuleTriangleTopology> triangles;
  std::vector<CapsuleTriangleReference> triRef;
  CapsuleMeshState state;

  int numNodes() const { return static_cast<int>(nodeTopo.size()); }
  int numEdges() const { return static_cast<int>(edges.size()); }
  int numTriangles() const { return static_cast<int>(triangles.size()); }
};

} // namespace Models
} // namespace ELFF
