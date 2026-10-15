#pragma once

#include "elff/models/capsule/CapsuleBending.hpp"
#include "elff/models/capsule/CapsuleElasticity.hpp"
#include "elff/models/capsule/CapsuleMeshBuilder.hpp"
#include "elff/models/capsule/CapsuleVolumeConstraint.hpp"

#include <memory>
#include <stdexcept>
#include <utility>

namespace ELFF {
namespace Models {

class Capsule {
public:
  explicit Capsule(CapsuleMesh mesh = {},
                   std::unique_ptr<CapsuleMembraneLaw> law =
                       std::make_unique<CapsuleNeoHookeanLaw>())
      : mesh_(std::move(mesh)), law_(std::move(law)) {
    if (!law_)
      throw std::invalid_argument("Capsule requires a membrane law");
  }

  Capsule(const Capsule &other)
      : mesh_(other.mesh_),
        law_(other.law_->clone()),
        bendingLaw_(other.bendingLaw_ ? other.bendingLaw_->clone() : nullptr),
        volumeConservationEnabled_(other.volumeConservationEnabled_) {}

  Capsule &operator=(const Capsule &other) {
    if (this == &other)
      return *this;
    mesh_ = other.mesh_;
    law_ = other.law_->clone();
    bendingLaw_ = other.bendingLaw_ ? other.bendingLaw_->clone() : nullptr;
    volumeConservationEnabled_ = other.volumeConservationEnabled_;
    return *this;
  }

  Capsule(Capsule &&) noexcept = default;
  Capsule &operator=(Capsule &&) noexcept = default;

  CapsuleMesh &mesh() { return mesh_; }
  const CapsuleMesh &mesh() const { return mesh_; }

  CapsuleMembraneLaw &law() { return *law_; }
  const CapsuleMembraneLaw &law() const { return *law_; }

  bool hasBendingLaw() const { return bendingLaw_ != nullptr; }

  CapsuleBendingLaw &bendingLaw() {
    if (!bendingLaw_)
      throw std::logic_error("Capsule does not have a bending law");
    return *bendingLaw_;
  }

  const CapsuleBendingLaw &bendingLaw() const {
    if (!bendingLaw_)
      throw std::logic_error("Capsule does not have a bending law");
    return *bendingLaw_;
  }

  void setLaw(std::unique_ptr<CapsuleMembraneLaw> law) {
    if (!law)
      throw std::invalid_argument("Capsule requires a membrane law");
    law_ = std::move(law);
  }

  void setBendingLaw(std::unique_ptr<CapsuleBendingLaw> law) {
    bendingLaw_ = std::move(law);
  }

  void initializeReferenceConfiguration() {
    CapsuleElasticityAssembler(law()).initializeReferenceConfiguration(mesh_);
    if (bendingLaw_)
      initializeReferenceCurvature();
  }

  void initializeReferenceCurvature() {
    if (!bendingLaw_)
      throw std::logic_error("Capsule requires a bending law");
    CapsuleBendingAssembler(bendingLaw()).initializeReferenceCurvature(mesh_);
  }

  void setConstantReferenceCurvature(double referenceCurvature) {
    if (!bendingLaw_)
      throw std::logic_error("Capsule requires a bending law");
    CapsuleBendingAssembler(bendingLaw()).setConstantReferenceCurvature(
        mesh_, referenceCurvature);
  }

  void setCurrentVolumeAsReference() {
    CapsuleVolumeConstraint::setCurrentVolumeAsReference(mesh_);
  }

  void initializeStructuralState(bool setVolumeReference = true) {
    initializeReferenceConfiguration();
    if (setVolumeReference)
      setCurrentVolumeAsReference();
  }

  void updateGeometry() {
    CapsuleGeometryOps::updateTriangleGeometry(mesh_);
    CapsuleGeometryOps::updateNodeNormals(mesh_);
    CapsuleGeometryOps::updateVolume(mesh_);
  }

  void computeElasticForces() {
    CapsuleElasticityAssembler(law()).assemble(mesh_);
    if (bendingLaw_)
      CapsuleBendingAssembler(bendingLaw()).assemble(mesh_);
  }

  double enforceVolumeConservation() {
    return CapsuleVolumeConstraint::enforce(mesh_);
  }

  void setVolumeConservationEnabled(bool enabled) {
    volumeConservationEnabled_ = enabled;
  }

  bool volumeConservationEnabled() const { return volumeConservationEnabled_; }

  void postAdvanceCorrection() {
    if (volumeConservationEnabled_)
      enforceVolumeConservation();
    else
      updateGeometry();
  }

private:
  CapsuleMesh mesh_;
  std::unique_ptr<CapsuleMembraneLaw> law_;
  std::unique_ptr<CapsuleBendingLaw> bendingLaw_;
  bool volumeConservationEnabled_ = true;
};

} // namespace Models
} // namespace ELFF
