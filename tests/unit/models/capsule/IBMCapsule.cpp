#include <gtest/gtest.h>

#include "elff/models/capsule/IBMCapsule.hpp"
#include "elff/models/capsule/CapsuleMeshBuilder.hpp"

#include <memory>
#include <vector>

namespace ELFF {
namespace Models {
namespace {

TEST(IBMCapsuleTest, MirrorsCapsuleMeshToCurrentIBMesh) {
  IBMCapsule capsule(CapsuleMeshBuilder::sphere({ 1.0, Vec3(1.0, 2.0, 3.0), 0 }));

  IBMesh &ibMesh = capsule.GetCurrent();

  ASSERT_EQ(ibMesh.GetNumberOfPoints(), capsule.mesh().numNodes());
  EXPECT_DOUBLE_EQ(ibMesh.GetPoints()[0].x, capsule.mesh().state.x(0, 0));
  EXPECT_DOUBLE_EQ(ibMesh.GetPoints()[0].y, capsule.mesh().state.x(1, 0));
  EXPECT_DOUBLE_EQ(ibMesh.GetPoints()[0].z, capsule.mesh().state.x(2, 0));
  EXPECT_GT(ibMesh.GetMeasures()[0], 0.0);
}

TEST(IBMCapsuleTest, MidpointComputesStructuralForces) {
  IBMCapsule capsule(CapsuleMeshBuilder::sphere({ 1.0, Vec3::Zero(), 1 }));

  std::vector<IBMesh::IBVertex> velocity(capsule.GetNumberOfPoints(),
                                         { 0.0, 0.0, 0.0 });
  velocity[0] = { 1.0, 0.0, 0.0 };

  IBMesh &midpoint = capsule.GetMidpoint(velocity, 0.1);

  EXPECT_TRUE(capsule.mesh().state.f.array().isFinite().all());
  EXPECT_GT(capsule.mesh().state.f.norm(), 0.0);
  const double measure = midpoint.GetMeasures()[0];
  EXPECT_GT(measure, 0.0);
  EXPECT_DOUBLE_EQ(midpoint.GetForces()[0].x,
                   capsule.mesh().state.f(0, 0) / measure);
  EXPECT_DOUBLE_EQ(midpoint.GetForces()[0].y,
                   capsule.mesh().state.f(1, 0) / measure);
  EXPECT_DOUBLE_EQ(midpoint.GetForces()[0].z,
                   capsule.mesh().state.f(2, 0) / measure);
}

TEST(IBMCapsuleTest, MembraneForceOpposesOutwardPerturbation) {
  IBMCapsule capsule(CapsuleMeshBuilder::sphere({ 1.0, Vec3::Zero(), 1 }));
  const Vec3 initialPosition = capsule.mesh().state.x.col(0);
  const Vec3 outward = initialPosition.normalized();

  capsule.mesh().state.x.col(0) += 0.1 * outward;
  capsule.computeElasticForces();

  EXPECT_LT(capsule.mesh().state.f.col(0).dot(outward), 0.0);
}

TEST(IBMCapsuleTest, NextAdvancesCurrentCapsulePositionFromVelocity) {
  IBMCapsule capsule(CapsuleMeshBuilder::sphere({ 1.0, Vec3::Zero(), 0 }));
  const Vec3 initialPosition = capsule.mesh().state.x.col(0);

  std::vector<IBMesh::IBVertex> velocity(capsule.GetNumberOfPoints(),
                                         { 1.0, 2.0, 3.0 });

  capsule.GetNext(velocity, 0.25);

  EXPECT_DOUBLE_EQ(capsule.mesh().state.x(0, 0), initialPosition.x() + 0.25);
  EXPECT_DOUBLE_EQ(capsule.mesh().state.x(1, 0), initialPosition.y() + 0.5);
  EXPECT_DOUBLE_EQ(capsule.mesh().state.x(2, 0), initialPosition.z() + 0.75);
  EXPECT_DOUBLE_EQ(capsule.mesh().state.v(0, 0), 1.0);
  EXPECT_DOUBLE_EQ(capsule.mesh().state.v(1, 0), 2.0);
  EXPECT_DOUBLE_EQ(capsule.mesh().state.v(2, 0), 3.0);
}

TEST(IBMCapsuleTest, PacksAndRestoresStructuralState) {
  IBMCapsule capsule(CapsuleMeshBuilder::sphere({ 1.0, Vec3::Zero(), 1 }));
  capsule.setLaw(std::make_unique<CapsuleSkalakLaw>(2.0, 15.0));
  capsule.setBendingLaw(std::make_unique<CapsuleHelfrichBendingLaw>(0.03));
  capsule.setConstantReferenceCurvature(-2.09);
  capsule.setVolumeConservationEnabled(true);

  std::vector<IBMesh::IBVertex> velocity(capsule.GetNumberOfPoints(),
                                         { 0.0, 0.0, 0.0 });
  velocity[0] = { 1.0, 2.0, 3.0 };
  capsule.GetNext(velocity, 0.25);
  capsule.computeElasticForces();

  IBModelState state;
  capsule.pack_state(state);
  ASSERT_EQ(state.ints.size(), 7u);
  EXPECT_EQ(state.ints[1], capsule.mesh().numNodes());
  EXPECT_GT(state.reals.size(), 0u);

  IBMCapsule restored(CapsuleMeshBuilder::sphere({ 1.0, Vec3::Zero(), 1 }));
  restored.unpack_state(state);

  EXPECT_TRUE(restored.volumeConservationEnabled());
  EXPECT_EQ(restored.mesh().numNodes(), capsule.mesh().numNodes());
  EXPECT_EQ(restored.mesh().numEdges(), capsule.mesh().numEdges());
  EXPECT_EQ(restored.mesh().numTriangles(), capsule.mesh().numTriangles());
  EXPECT_TRUE(restored.mesh().state.x.isApprox(capsule.mesh().state.x));
  EXPECT_TRUE(restored.mesh().state.v.isApprox(capsule.mesh().state.v));
  EXPECT_TRUE(restored.mesh().state.f.isApprox(capsule.mesh().state.f));
  EXPECT_TRUE(restored.mesh().state.refCurv.isApprox(capsule.mesh().state.refCurv));
  EXPECT_DOUBLE_EQ(restored.mesh().state.initialVolume,
                   capsule.mesh().state.initialVolume);

  auto *law = dynamic_cast<CapsuleSkalakLaw *>(&restored.law());
  ASSERT_NE(law, nullptr);
  EXPECT_DOUBLE_EQ(law->Es, 2.0);
  EXPECT_DOUBLE_EQ(law->C, 15.0);
  auto *bending =
      dynamic_cast<CapsuleHelfrichBendingLaw *>(&restored.bendingLaw());
  ASSERT_NE(bending, nullptr);
  EXPECT_DOUBLE_EQ(bending->Eb, 0.03);

  IBMesh &ibMesh = restored.GetCurrent();
  EXPECT_DOUBLE_EQ(ibMesh.GetPoints()[0].x, restored.mesh().state.x(0, 0));
  EXPECT_DOUBLE_EQ(ibMesh.GetVelocity()[0].y, restored.mesh().state.v(1, 0));
  EXPECT_GT(ibMesh.GetMeasures()[0], 0.0);
  EXPECT_DOUBLE_EQ(ibMesh.GetForces()[0].x,
                   restored.mesh().state.f(0, 0) / ibMesh.GetMeasures()[0]);
}

} // namespace
} // namespace Models
} // namespace ELFF
