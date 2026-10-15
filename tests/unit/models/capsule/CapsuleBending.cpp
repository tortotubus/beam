#include <gtest/gtest.h>

#include "elff/models/capsule/Capsule.hpp"
#include "elff/models/capsule/CapsuleBending.hpp"
#include "elff/models/capsule/CapsuleMeshBuilder.hpp"

#include <memory>

namespace ELFF {
namespace Models {
namespace {

TEST(CapsuleBendingLawTest, LinearLawUsesOnlyCurvatureLaplacian) {
  CapsuleLinearBendingLaw law(0.25);

  EXPECT_DOUBLE_EQ(law.surfaceForceDensity(1.0, 0.5, 0.25, 3.0), 1.5);
}

TEST(CapsuleBendingLawTest, HelfrichLawIncludesNonlinearCurvatureTerms) {
  CapsuleHelfrichBendingLaw law(0.25);

  const double density = law.surfaceForceDensity(2.0, 0.5, 0.25, 3.0);
  const double expected =
      2.0 * 0.25 *
      (3.0 + 2.0 * (2.0 - 0.5) * (2.0 * 2.0 - 0.25 + 0.5 * 2.0));

  EXPECT_DOUBLE_EQ(density, expected);
}

TEST(CapsuleBendingAssemblerTest,
     ReferenceCurvatureMakesInitialConfigurationStressFree) {
  CapsuleMesh mesh = CapsuleMeshBuilder::sphere({ 1.0, Vec3::Zero(), 1 });
  CapsuleLinearBendingLaw law(0.1);
  CapsuleBendingAssembler assembler(law);

  assembler.initializeReferenceCurvature(mesh);
  mesh.state.zeroForces();
  assembler.assemble(mesh);

  EXPECT_NEAR(mesh.state.f.norm(), 0.0, 1.0e-12);
}

TEST(CapsuleBendingAssemblerTest, PerturbedMeshProducesFiniteBendingForces) {
  CapsuleMesh mesh = CapsuleMeshBuilder::sphere({ 1.0, Vec3::Zero(), 1 });
  CapsuleHelfrichBendingLaw law(0.1);
  CapsuleBendingAssembler assembler(law);

  assembler.initializeReferenceCurvature(mesh);
  mesh.state.x.col(0) *= 1.1;
  mesh.state.zeroForces();
  assembler.assemble(mesh);

  EXPECT_TRUE(mesh.state.f.array().isFinite().all());
  EXPECT_GT(mesh.state.f.norm(), 0.0);
}

TEST(CapsuleBendingAssemblerTest, MixedNodeAreasRemainPositiveAndSurfaceScaled) {
  for (CapsuleMesh mesh :
       { CapsuleMeshBuilder::sphere({ 1.0, Vec3::Zero(), 2 }),
         CapsuleMeshBuilder::biconcave({ 1.0, Vec3::Zero(), 2 }) }) {
    CapsuleGeometryOps::updateTriangleGeometry(mesh);

    const std::vector<double> nodeArea =
        CapsuleBendingAssembler::computeNodeAreas(mesh);

    double nodeAreaSum = 0.0;
    for (const double area : nodeArea) {
      EXPECT_GT(area, 0.0);
      nodeAreaSum += area;
    }

    double triangleAreaSum = 0.0;
    for (const auto &geom : mesh.state.triGeom)
      triangleAreaSum += geom.area;

    EXPECT_NEAR(nodeAreaSum, triangleAreaSum, 0.01 * triangleAreaSum);
  }
}

TEST(CapsuleBendingAssemblerTest, CapsuleAddsOptionalBendingForces) {
  Capsule capsule(CapsuleMeshBuilder::sphere({ 1.0, Vec3::Zero(), 1 }));
  capsule.setBendingLaw(std::make_unique<CapsuleLinearBendingLaw>(0.1));
  capsule.initializeStructuralState();

  Capsule referenceCapsule = capsule;
  referenceCapsule.computeElasticForces();
  EXPECT_NEAR(referenceCapsule.mesh().state.f.norm(), 0.0, 1.0e-12);

  capsule.mesh().state.x.col(0) *= 1.1;
  capsule.computeElasticForces();
  EXPECT_TRUE(capsule.mesh().state.f.array().isFinite().all());
  EXPECT_GT(capsule.mesh().state.f.norm(), 0.0);
}

} // namespace
} // namespace Models
} // namespace ELFF
