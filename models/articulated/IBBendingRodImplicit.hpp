#pragma once

#include "elff/models/articulated/BendingRodImplicit.hpp"
#include "elff/models/ibm/IBForceCoupled.hpp"

namespace ELFF::Models {
/** Immersed-boundary adapter mirroring the common rod mesh onto IB points. */
class IBBendingRodImplicit : public IBForceCoupled, public BendingRodImplicit {
public:
  explicit IBBendingRodImplicit(const BendingRodImplicitConfig &config);
  const BendingRodImplicit &solver() const { return *this; }
  using BendingRodImplicit::solve;
  void solve(real_t dt, const std::vector<Vector>& density) override;
  void ComputeNextPoints(std::vector<IBMesh::IBVertex> force, real_t dt) override;
  void pack_state(IBModelState &s) const override;
  void unpack_state(const IBModelState &s) override;

protected:
  void setState(const Eigen::VectorXd &state) override;
  void syncMeshes();
};
} // namespace ELFF::Models
