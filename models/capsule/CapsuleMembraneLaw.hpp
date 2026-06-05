#pragma once

#include "elff/models/capsule/CapsuleMesh.hpp"

#include <cmath>
#include <memory>

namespace ELFF {
namespace Models {

class CapsuleMembraneLaw {
public:
  virtual ~CapsuleMembraneLaw() = default;
  virtual Vec2 dWdLambda(const Vec2& lambda) const = 0;
  virtual std::unique_ptr<CapsuleMembraneLaw> clone() const = 0;
};

class CapsuleNeoHookeanLaw final : public CapsuleMembraneLaw {
public:
  explicit CapsuleNeoHookeanLaw(double elasticModulus = 1.0)
    : Es(elasticModulus) {}

  double Es = 1.0;

  Vec2 dWdLambda(const Vec2& lambda) const override {
    const double l1 = lambda[0], l2 = lambda[1];
    return {
      Es / (3.0 * l1) * (l1 * l1 - 1.0 / std::pow(l1 * l2, 2)),
      Es / (3.0 * l2) * (l2 * l2 - 1.0 / std::pow(l1 * l2, 2))
    };
  }

  std::unique_ptr<CapsuleMembraneLaw> clone() const override {
    return std::make_unique<CapsuleNeoHookeanLaw>(*this);
  }
};

class CapsuleSkalakLaw final : public CapsuleMembraneLaw {
public:
  CapsuleSkalakLaw(double elasticModulus = 1.0,
                   double areaDilatationModulus = 1.0)
    : Es(elasticModulus), C(areaDilatationModulus) {}

  double Es = 1.0;
  double C = 1.0;

  Vec2 dWdLambda(const Vec2& lambda) const override {
    const double l1 = lambda[0], l2 = lambda[1];
    const double J2m1 = std::pow(l1 * l2, 2) - 1.0;
    return {
      Es * (l1 * (l1 * l1 - 1.0) + C * l1 * l2 * l2 * J2m1),
      Es * (l2 * (l2 * l2 - 1.0) + C * l2 * l1 * l1 * J2m1)
    };
  }

  std::unique_ptr<CapsuleMembraneLaw> clone() const override {
    return std::make_unique<CapsuleSkalakLaw>(*this);
  }
};

} // namespace Models
} // namespace ELFF
