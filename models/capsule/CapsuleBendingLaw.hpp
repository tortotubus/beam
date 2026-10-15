#pragma once

#include "elff/models/capsule/CapsuleMesh.hpp"

#include <memory>

namespace ELFF {
namespace Models {

class CapsuleBendingLaw {
public:
  virtual ~CapsuleBendingLaw() = default;

  virtual double surfaceForceDensity(double meanCurvature,
                                     double referenceCurvature,
                                     double gaussianCurvature,
                                     double curvatureLaplacian) const = 0;

  virtual std::unique_ptr<CapsuleBendingLaw> clone() const = 0;
};

class CapsuleLinearBendingLaw final : public CapsuleBendingLaw {
public:
  explicit CapsuleLinearBendingLaw(double bendingModulus = 0.0)
      : Eb(bendingModulus) {}

  double Eb = 0.0;

  double surfaceForceDensity(double,
                             double,
                             double,
                             double curvatureLaplacian) const override {
    return 2.0 * Eb * curvatureLaplacian;
  }

  std::unique_ptr<CapsuleBendingLaw> clone() const override {
    return std::make_unique<CapsuleLinearBendingLaw>(*this);
  }
};

class CapsuleHelfrichBendingLaw final : public CapsuleBendingLaw {
public:
  explicit CapsuleHelfrichBendingLaw(double bendingModulus = 0.0)
      : Eb(bendingModulus) {}

  double Eb = 0.0;

  double surfaceForceDensity(double meanCurvature,
                             double referenceCurvature,
                             double gaussianCurvature,
                             double curvatureLaplacian) const override {
    const double curvatureDifference = meanCurvature - referenceCurvature;
    return 2.0 * Eb *
           (curvatureLaplacian +
            2.0 * curvatureDifference *
                (meanCurvature * meanCurvature - gaussianCurvature +
                 referenceCurvature * meanCurvature));
  }

  std::unique_ptr<CapsuleBendingLaw> clone() const override {
    return std::make_unique<CapsuleHelfrichBendingLaw>(*this);
  }
};

} // namespace Models
} // namespace ELFF
