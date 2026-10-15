#pragma once

#include "elff/models/beam/EulerBeamInextensibleADDM.hpp"

namespace ELFF {
namespace Models {

/**
 * @brief Thin test harness for experiments based on EulerBeamInextensibleADDM.
 *
 * This class intentionally inherits the parent's constructors and behavior
 * unchanged. It gives benchmark/reference-comparison work a distinct type that
 * can diverge from the production ADDM implementation later without changing
 * existing ADDM call sites.
 */
class EulerBeamInextensibleADDMTestHarness
    : public EulerBeamInextensibleADDM {
public:
  using EulerBeamInextensibleADDM::EulerBeamInextensibleADDM;

  /**
   * @brief Advance one average-acceleration Newmark step for the unconstrained
   * linear beam equation.
   *
   * The supplied load is the Newmark averaged load over the interval, i.e.
   * f_tilde = 0.5 * (f^n + f^{n+1}) for the manufactured-solution benchmark.
   * This intentionally bypasses the ADDM projection, multiplier update, and
   * penalty terms so the test harness can reproduce the linear case in
   * Basting et al.
  */
  void solve_linear(real_t dt,
                    const std::vector<std::array<real_t, 3>> &averaged_load);
};

} // namespace Models
} // namespace ELFF
