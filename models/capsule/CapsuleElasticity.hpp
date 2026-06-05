#pragma once

#include "elff/models/capsule/CapsuleGeometry.hpp"
#include "elff/models/capsule/CapsuleMembraneLaw.hpp"

#include <cmath>
#include <stdexcept>

namespace ELFF {
namespace Models {

class CapsuleElasticityAssembler {
public:
  explicit CapsuleElasticityAssembler(const CapsuleMembraneLaw& law)
    : law_(law) {}

  void initializeReferenceConfiguration(CapsuleMesh& mesh) const {
    CapsuleGeometryOps::updateTriangleGeometry(mesh);
    mesh.triRef.resize(mesh.numTriangles());

    for (int tid = 0; tid < mesh.numTriangles(); ++tid) {
      const auto refPlane =
        CapsuleGeometryOps::rotateTriangleToReferencePlane(mesh, tid);
      const Mat32& refNodes = refPlane.first;

      auto& ref = mesh.triRef[tid];
      ref.refShape.col(0) = refNodes.col(0).head<2>();
      ref.refShape.col(1) = refNodes.col(1).head<2>();
      ref.restArea = mesh.state.triGeom[tid].area;

      const double x1 = ref.refShape(0,0), y1 = ref.refShape(1,0);
      const double x2 = ref.refShape(0,1), y2 = ref.refShape(1,1);
      const double det = x1 * y2 - y1 * x2;
      if (std::abs(det) <= 1.0e-14)
        throw std::runtime_error("Capsule reference triangle is degenerate");

      ref.dNdXi.row(0) << (y1 - y2) / det, (x2 - x1) / det;
      ref.dNdXi.row(1) << y2 / det, -x2 / det;
      ref.dNdXi.row(2) << -y1 / det, x1 / det;
    }
  }

  void assemble(CapsuleMesh& mesh) const {
    CapsuleGeometryOps::updateTriangleGeometry(mesh);
    CapsuleGeometryOps::updateNodeNormals(mesh);
    mesh.state.zeroForces();

    for (int tid = 0; tid < mesh.numTriangles(); ++tid) {
      const auto& tri = mesh.triangles[tid];
      const auto& ref = mesh.triRef[tid];
      auto& geom = mesh.state.triGeom[tid];

      const auto [curNodesInRefPlane, RrefToCur] =
        CapsuleGeometryOps::rotateTriangleToReferencePlane(mesh, tid);

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
      if (lambda[0] <= 0.0 || lambda[1] <= 0.0)
        throw std::runtime_error(
          "Capsule triangle has a degenerate deformation gradient");

      Vec2 dW = law_.dWdLambda(lambda);
      geom.stretch = lambda;
      geom.tension << dW[0] / lambda[1], dW[1] / lambda[0];

      for (int a = 0; a < 3; ++a) {
        Vec2 fPlane = localNodalForce(a, F, C, lambda, ref.dNdXi, dW);
        Vec3 f3 = geom.area * RrefToCur.leftCols<2>() * fPlane;
        mesh.state.f.col(tri.nodes[a]) -= f3;
      }
    }
  }

private:
  const CapsuleMembraneLaw& law_;

  static Vec2 localNodalForce(
    int a,
    const Mat2& F,
    const Mat2& C,
    const Vec2& lambda,
    const Eigen::Matrix<double,3,2>& dNdXi,
    const Vec2& dW) {
    const double den = std::sqrt(std::pow(C(0,0) - C(1,1), 2) +
                                 4.0 * std::pow(C(0,1), 2));
    const int sign[2] = {-1, 1};
    const Vec2 sf = dNdXi.row(a).transpose();

    Eigen::Matrix<double,2,2> dLambdaDv =
      Eigen::Matrix<double,2,2>::Zero();
    for (int k = 0; k < 2; ++k) {
      if (lambda[k] <= 0.0)
        continue;

      for (int l = 0; l < 2; ++l) {
        const int c1 = l;
        const int c2 = (l + 1) % 2;

        dLambdaDv(k,l) =
          (sf[c1] * F(c1,c1) + sf[c2] * F(c1,c2)) /
          (2.0 * lambda[k]);

        if (den > 1.0e-12) {
          dLambdaDv(k,l) += sign[k] *
            ((C(c1,c1) - C(c2,c2)) *
               (sf[c1] * F(c1,c1) - sf[c2] * F(c1,c2)) +
             2.0 * C(0,1) *
               (sf[c1] * F(c1,c2) + sf[c2] * F(c1,c1))) /
            (den * 2.0 * lambda[k]);
        }
      }
    }

    return {
      dW[0] * dLambdaDv(0,0) + dW[1] * dLambdaDv(1,0),
      dW[0] * dLambdaDv(0,1) + dW[1] * dLambdaDv(1,1)
    };
  }
};

} // namespace Models
} // namespace ELFF
