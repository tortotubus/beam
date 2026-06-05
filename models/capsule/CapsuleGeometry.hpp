#pragma once

#include "elff/models/capsule/CapsuleMesh.hpp"

#include <stdexcept>
#include <utility>

namespace ELFF {
namespace Models {

class CapsuleGeometryOps {
public:
  static void updateTriangleGeometry(CapsuleMesh& mesh) {
    if (mesh.state.triGeom.size() != static_cast<size_t>(mesh.numTriangles()))
      mesh.state.triGeom.resize(mesh.numTriangles());

    for (int tid = 0; tid < mesh.numTriangles(); ++tid) {
      const auto& tri = mesh.triangles[tid];
      const Vec3 x0 = mesh.state.x.col(tri.nodes[0]);
      const Vec3 d10 = mesh.state.x.col(tri.nodes[1]) - x0;
      const Vec3 d20 = mesh.state.x.col(tri.nodes[2]) - x0;
      const Vec3 areaNormal = d10.cross(d20);
      const double norm = areaNormal.norm();

      auto& geom = mesh.state.triGeom[tid];
      geom.area = 0.5 * norm;
      geom.normal = Vec3::Zero();
      if (norm > 0.0)
        geom.normal = areaNormal / norm;
      geom.centroid = x0 + (d10 + d20) / 3.0;
    }
  }

  static void updateNodeNormals(CapsuleMesh& mesh) {
    mesh.state.normal.setZero(3, mesh.numNodes());
    for (int tid = 0; tid < mesh.numTriangles(); ++tid) {
      const auto& tri = mesh.triangles[tid];
      const auto& geom = mesh.state.triGeom[tid];
      for (int a = 0; a < 3; ++a)
        mesh.state.normal.col(tri.nodes[a]) += geom.area * geom.normal;
    }

    for (int nid = 0; nid < mesh.numNodes(); ++nid) {
      const double norm = mesh.state.normal.col(nid).norm();
      if (norm > 0.0)
        mesh.state.normal.col(nid) /= norm;
    }
  }

  static void updateCentroid(CapsuleMesh& mesh) {
    if (mesh.numNodes() == 0) {
      mesh.state.centroid = Vec3::Zero();
      return;
    }

    const Vec3 x0 = mesh.state.x.col(0);
    Vec3 sum = Vec3::Zero();
    for (int nid = 0; nid < mesh.numNodes(); ++nid)
      sum += mesh.state.x.col(nid) - x0;
    mesh.state.centroid = x0 + sum / static_cast<double>(mesh.numNodes());
  }

  static void updateVolume(CapsuleMesh& mesh) {
    updateCentroid(mesh);
    double volume = 0.0;

    for (int tid = 0; tid < mesh.numTriangles(); ++tid) {
      const auto& tri = mesh.triangles[tid];
      const Vec3 r0 = mesh.state.x.col(tri.nodes[0]) - mesh.state.centroid;
      const Vec3 r1 = mesh.state.x.col(tri.nodes[1]) - mesh.state.centroid;
      const Vec3 r2 = mesh.state.x.col(tri.nodes[2]) - mesh.state.centroid;
      volume += r0.dot(r1.cross(r2)) / 6.0;
    }

    mesh.state.volume = volume;
  }

  static double edgeLength(const CapsuleMesh& mesh, int eid) {
    const auto& edge = mesh.edges[eid];
    return (mesh.state.x.col(edge.nodes[0]) -
            mesh.state.x.col(edge.nodes[1])).norm();
  }

  static std::pair<Mat32, Mat3> rotateTriangleToReferencePlane(
    const CapsuleMesh& mesh, int tid) {
    const auto& tri = mesh.triangles[tid];
    const auto& geom = mesh.state.triGeom[tid];

    const Vec3 x0 = mesh.state.x.col(tri.nodes[0]);
    const Vec3 d10 = mesh.state.x.col(tri.nodes[1]) - x0;
    const Vec3 d20 = mesh.state.x.col(tri.nodes[2]) - x0;

    const double e0Norm = d20.norm();
    if (e0Norm <= 0.0)
      throw std::runtime_error(
        "Capsule triangle has a degenerate node0-node2 edge");

    Mat3 currentBasis = Mat3::Zero();
    currentBasis.col(0) = d20 / e0Norm;
    currentBasis.col(2) = -geom.normal;
    currentBasis.col(1) = currentBasis.col(2).cross(currentBasis.col(0));

    const double e1Norm = currentBasis.col(1).norm();
    if (e1Norm <= 0.0)
      throw std::runtime_error("Capsule triangle has a degenerate basis");
    currentBasis.col(1) /= e1Norm;

    const Mat3 curToRef = currentBasis.transpose();
    const Mat3 refToCur = currentBasis;

    Mat32 local = Mat32::Zero();
    local.col(0) = curToRef * d10;
    local.col(1) = curToRef * d20;
    return {local, refToCur};
  }
};

} // namespace Models
} // namespace ELFF
