#pragma once

#include "elff/models/capsule/CapsuleBendingLaw.hpp"
#include "elff/models/capsule/CapsuleGeometry.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace ELFF {
namespace Models {

class CapsuleBendingAssembler {
public:
  explicit CapsuleBendingAssembler(const CapsuleBendingLaw &law)
      : law_(law) {}

  void initializeReferenceCurvature(CapsuleMesh &mesh) const {
    updateCurvatures(mesh);
    mesh.state.refCurv = mesh.state.meanCurv;
  }

  void setConstantReferenceCurvature(CapsuleMesh &mesh,
                                     double referenceCurvature) const {
    ensureNodeState(mesh);
    mesh.state.refCurv.setConstant(referenceCurvature);
  }

  void assemble(CapsuleMesh &mesh) const {
    updateCurvatures(mesh);

    const std::vector<double> nodeArea = computeNodeAreas(mesh);
    const Eigen::VectorXd curvatureLaplacian =
        computeCurvatureLaplacian(mesh, nodeArea);

    for (int nid = 0; nid < mesh.numNodes(); ++nid) {
      const double traction = law_.surfaceForceDensity(
          mesh.state.meanCurv(nid),
          mesh.state.refCurv(nid),
          mesh.state.gaussCurv(nid),
          curvatureLaplacian(nid));
      mesh.state.f.col(nid) +=
          nodeArea[static_cast<size_t>(nid)] * traction *
          mesh.state.normal.col(nid);
    }
  }

  static void updateCurvatures(CapsuleMesh &mesh) {
    constexpr double pi = 3.141592653589793238462643383279502884;

    CapsuleGeometryOps::updateTriangleGeometry(mesh);
    CapsuleGeometryOps::updateNodeNormals(mesh);
    ensureNodeState(mesh);

    const std::vector<double> nodeArea = computeNodeAreas(mesh);
    const std::vector<double> angleSum = computeAngleSums(mesh);
    const Eigen::Matrix3Xd meanCurvatureNormal =
        computeMeanCurvatureNormal(mesh);

    for (int nid = 0; nid < mesh.numNodes(); ++nid) {
      const double area = nodeArea[static_cast<size_t>(nid)];
      if (area <= 0.0)
        throw std::runtime_error("Capsule node has zero bending area");

      const Vec3 normal = mesh.state.normal.col(nid);
      mesh.state.meanCurv(nid) =
          -0.5 * meanCurvatureNormal.col(nid).dot(normal);
      mesh.state.gaussCurv(nid) =
          (2.0 * pi - angleSum[static_cast<size_t>(nid)]) / area;
    }
  }

  static std::vector<double> computeNodeAreas(const CapsuleMesh &mesh) {
    std::vector<double> nodeArea(static_cast<size_t>(mesh.numNodes()), 0.0);

    for (int tid = 0; tid < mesh.numTriangles(); ++tid) {
      const auto &tri = mesh.triangles[tid];
      const double area = mesh.state.triGeom[tid].area / 3.0;
      for (const int nid : tri.nodes)
        nodeArea[static_cast<size_t>(nid)] += area;
    }

    return nodeArea;
  }

private:
  const CapsuleBendingLaw &law_;

  static void ensureNodeState(CapsuleMesh &mesh) {
    const int nNodes = mesh.numNodes();
    if (mesh.state.normal.cols() != nNodes)
      mesh.state.normal.setZero(3, nNodes);
    if (mesh.state.meanCurv.size() != nNodes)
      mesh.state.meanCurv.setZero(nNodes);
    if (mesh.state.refCurv.size() != nNodes)
      mesh.state.refCurv.setZero(nNodes);
    if (mesh.state.gaussCurv.size() != nNodes)
      mesh.state.gaussCurv.setZero(nNodes);
  }

  static double cotangent(const Vec3 &a, const Vec3 &b) {
    const double crossNorm = a.cross(b).norm();
    if (crossNorm <= 1.0e-14)
      return 0.0;
    return a.dot(b) / crossNorm;
  }

  static std::vector<double> computeAngleSums(const CapsuleMesh &mesh) {
    std::vector<double> angleSum(static_cast<size_t>(mesh.numNodes()), 0.0);

    for (const auto &tri : mesh.triangles) {
      for (int a = 0; a < 3; ++a) {
        const int nid = tri.nodes[a];
        const Vec3 x0 = mesh.state.x.col(nid);
        const Vec3 e1 = mesh.state.x.col(tri.nodes[(a + 1) % 3]) - x0;
        const Vec3 e2 = mesh.state.x.col(tri.nodes[(a + 2) % 3]) - x0;
        const double denom = e1.norm() * e2.norm();
        if (denom <= 0.0)
          continue;
        const double cosine =
            std::clamp(e1.dot(e2) / denom, -1.0, 1.0);
        angleSum[static_cast<size_t>(nid)] += std::acos(cosine);
      }
    }

    return angleSum;
  }

  static Eigen::Matrix3Xd computeMeanCurvatureNormal(const CapsuleMesh &mesh) {
    Eigen::Matrix3Xd meanCurvatureNormal =
        Eigen::Matrix3Xd::Zero(3, mesh.numNodes());
    const std::vector<double> nodeArea = computeNodeAreas(mesh);

    for (const auto &edge : mesh.edges) {
      if (edge.triangles[0] < 0 || edge.triangles[1] < 0)
        continue;

      const int i = edge.nodes[0];
      const int j = edge.nodes[1];
      const double weight = cotangentOppositeEdge(mesh, edge, 0) +
                            cotangentOppositeEdge(mesh, edge, 1);
      const Vec3 edgeVector = mesh.state.x.col(j) - mesh.state.x.col(i);

      meanCurvatureNormal.col(i) += weight * edgeVector;
      meanCurvatureNormal.col(j) -= weight * edgeVector;
    }

    for (int nid = 0; nid < mesh.numNodes(); ++nid) {
      const double area = nodeArea[static_cast<size_t>(nid)];
      if (area > 0.0)
        meanCurvatureNormal.col(nid) /= 2.0 * area;
    }

    return meanCurvatureNormal;
  }

  static double cotangentOppositeEdge(const CapsuleMesh &mesh,
                                      const CapsuleEdge &edge,
                                      int incidentTriangleIndex) {
    const int tid = edge.triangles[incidentTriangleIndex];
    if (tid < 0)
      return 0.0;

    const auto &tri = mesh.triangles[tid];
    int opposite = -1;
    for (const int nid : tri.nodes) {
      if (nid != edge.nodes[0] && nid != edge.nodes[1]) {
        opposite = nid;
        break;
      }
    }

    if (opposite < 0)
      return 0.0;

    const Vec3 xOpposite = mesh.state.x.col(opposite);
    return cotangent(mesh.state.x.col(edge.nodes[0]) - xOpposite,
                     mesh.state.x.col(edge.nodes[1]) - xOpposite);
  }

  static Eigen::VectorXd computeCurvatureLaplacian(
      const CapsuleMesh &mesh,
      const std::vector<double> &nodeArea) {
    Eigen::VectorXd laplacian = Eigen::VectorXd::Zero(mesh.numNodes());

    for (const auto &edge : mesh.edges) {
      if (edge.triangles[0] < 0 || edge.triangles[1] < 0)
        continue;

      const int i = edge.nodes[0];
      const int j = edge.nodes[1];
      const double weight = 0.5 * (cotangentOppositeEdge(mesh, edge, 0) +
                                   cotangentOppositeEdge(mesh, edge, 1));
      const double diff =
          (mesh.state.meanCurv(j) - mesh.state.refCurv(j)) -
          (mesh.state.meanCurv(i) - mesh.state.refCurv(i));

      laplacian(i) += weight * diff;
      laplacian(j) -= weight * diff;
    }

    for (int nid = 0; nid < mesh.numNodes(); ++nid) {
      const double area = nodeArea[static_cast<size_t>(nid)];
      if (area > 0.0)
        laplacian(nid) /= area;
    }

    return laplacian;
  }
};

} // namespace Models
} // namespace ELFF
