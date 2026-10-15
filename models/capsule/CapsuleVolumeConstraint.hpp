#pragma once

#include "elff/models/capsule/CapsuleGeometry.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

namespace ELFF {
namespace Models {

class CapsuleVolumeConstraint {
public:
  static void setCurrentVolumeAsReference(CapsuleMesh& mesh) {
    CapsuleGeometryOps::updateVolume(mesh);
    mesh.state.initialVolume = mesh.state.volume;
  }

  static Vec3 computeAlpha(const CapsuleMesh& mesh, int nid) {
    Vec3 alpha = Vec3::Zero();
    const auto incidentTriangles = trianglesForNode(mesh, nid);

    for (const int tid : incidentTriangles) {
      const auto& tri = mesh.triangles[tid];
      int localId = 0;
      while (localId < 3 && tri.nodes[localId] != nid)
        ++localId;
      if (localId == 3)
        continue;

      const int n0 = tri.nodes[(localId + 1) % 3];
      const int n1 = tri.nodes[(localId + 2) % 3];
      const Vec3 x0 = position(mesh, n0);
      const Vec3 x1 = position(mesh, n1);
      alpha -= x0.cross(x1) / 12.0;
    }

    return alpha;
  }

  static double enforce(CapsuleMesh& mesh) {
    CapsuleGeometryOps::updateVolume(mesh);
    if (mesh.numNodes() == 0 || mesh.numTriangles() == 0)
      return 0.0;

    Eigen::Matrix3Xd alpha(3, mesh.numNodes());
    for (int nid = 0; nid < mesh.numNodes(); ++nid)
      alpha.col(nid) = computeAlpha(mesh, nid);

    std::array<double,4> coeff{0.0, 0.0, 0.0, 0.0};
    for (int tid = 0; tid < mesh.numTriangles(); ++tid) {
      const auto& tri = mesh.triangles[tid];
      for (int k = 0; k < 3; ++k) {
        const int i0 = tri.nodes[k];
        const int i1 = tri.nodes[(k + 1) % 3];
        const int i2 = tri.nodes[(k + 2) % 3];

        const Vec3 x0 = position(mesh, i0);
        const Vec3 x1 = position(mesh, i1);
        const Vec3 x2 = position(mesh, i2);
        const Vec3 a0 = alpha.col(i0);
        const Vec3 a1 = alpha.col(i1);
        const Vec3 a2 = alpha.col(i2);

        coeff[3] += a0.dot(a1.cross(a2));
        coeff[2] += x0.dot(a1.cross(a2));
        coeff[1] += a0.dot(x1.cross(x2));
      }
    }

    coeff[3] /= 18.0;
    coeff[2] /= 6.0;
    coeff[1] /= 6.0;
    coeff[0] = mesh.state.volume - mesh.state.initialVolume;

    double normalize = 0.0;
    for (const double c : coeff)
      normalize = std::max(normalize, std::abs(c));

    if (normalize <= 1.0e-30)
      return 0.0;

    for (double& c : coeff)
      c /= normalize;

    const double lambda = smallestRealRoot(coeff);
    for (int nid = 0; nid < mesh.numNodes(); ++nid)
      mesh.state.x.col(nid) += lambda * alpha.col(nid);

    CapsuleGeometryOps::updateTriangleGeometry(mesh);
    CapsuleGeometryOps::updateVolume(mesh);
    return lambda;
  }

private:
  static Vec3 position(const CapsuleMesh& mesh, int nid) {
    return mesh.state.x.col(nid);
  }

  static std::vector<int> trianglesForNode(const CapsuleMesh& mesh, int nid) {
    if (nid < static_cast<int>(mesh.nodeTopo.size()) &&
        !mesh.nodeTopo[nid].triangles.empty())
      return mesh.nodeTopo[nid].triangles;

    std::vector<int> triangles;
    for (int tid = 0; tid < mesh.numTriangles(); ++tid) {
      const auto& tri = mesh.triangles[tid];
      if (tri.nodes[0] == nid || tri.nodes[1] == nid || tri.nodes[2] == nid)
        triangles.push_back(tid);
    }
    return triangles;
  }

  static double smallestRealRoot(const std::array<double,4>& c) {
    constexpr double eps = 1.0e-10;
    if (std::abs(c[3]) > eps)
      return cubicSmallestRealRoot(c);
    if (std::abs(c[2]) > eps)
      return quadraticSmallestRealRoot(c[2], c[1], c[0]);
    if (std::abs(c[1]) > eps)
      return -c[0] / c[1];
    return 0.0;
  }

  static double quadraticSmallestRealRoot(double a, double b, double c) {
    constexpr double eps = 1.0e-16;
    if (std::abs(a) < eps)
      return std::abs(b) > eps ? -c / b : 0.0;

    const double delta = b * b - 4.0 * a * c;
    if (delta < -eps)
      return 0.0;
    if (std::abs(delta) <= eps)
      return -b / (2.0 * a);

    const double rootDelta = std::sqrt(delta);
    const double r0 = (-b - rootDelta) / (2.0 * a);
    const double r1 = (-b + rootDelta) / (2.0 * a);
    return std::abs(r0) < std::abs(r1) ? r0 : r1;
  }

  static double cubicSmallestRealRoot(const std::array<double,4>& c) {
    constexpr double eps = 1.0e-10;
    constexpr double pi = 3.141592653589793238462643383279502884;

    const double alpha = c[2] / c[3];
    const double beta = c[1] / c[3];
    const double gamma = c[0] / c[3];
    const double Q = (alpha * alpha - 3.0 * beta) / 9.0;
    const double R = (2.0 * alpha * alpha * alpha -
                      9.0 * alpha * beta + 27.0 * gamma) / 54.0;
    const double Q3 = Q * Q * Q;
    const double R2 = R * R;

    if (R2 < Q3) {
      const double theta = std::acos(R / std::sqrt(Q3));
      const double r0 =
        -2.0 * std::sqrt(Q) * std::cos(theta / 3.0) - alpha / 3.0;
      const double r1 =
        -2.0 * std::sqrt(Q) * std::cos((theta + 2.0 * pi) / 3.0) -
        alpha / 3.0;
      const double r2 =
        -2.0 * std::sqrt(Q) * std::cos((theta - 2.0 * pi) / 3.0) -
        alpha / 3.0;

      double result = r0;
      if (std::abs(r1) < std::abs(result))
        result = r1;
      if (std::abs(r2) < std::abs(result))
        result = r2;
      return result;
    }

    const double radical = std::sqrt(std::max(0.0, R2 - Q3));
    const double A = -std::copysign(std::cbrt(std::abs(R) + radical), R);
    const double B = std::abs(A) > eps ? Q / A : 0.0;
    return A + B - alpha / 3.0;
  }
};

} // namespace Models
} // namespace ELFF
