#include "elff/models/articulated/BendingRodMesh.hpp"
#include <cmath>
#include "elff/general/error.hpp"

namespace ELFF::Models {
BendingRodMesh::BendingRodMesh(size_t nodes, real_t length) : length_(length)
{
  if (nodes < 2 || !std::isfinite(length) || length <= 0)
    ELFF_ABORT("rod mesh requires at least two nodes and positive finite length");
  ds_ = length / (nodes - 1);
  if (!std::isfinite(ds_) || ds_ <= 0)
    ELFF_ABORT("invalid rod mesh spacing");
  s_.resize(nodes);
  centerline_.resize(nodes);
  velocity_.resize(nodes, Vector{0, 0, 0});
  for (size_t i = 0; i < nodes; ++i) {
    s_[i] = i == nodes - 1 ? length : i * ds_;
    centerline_[i] = {s_[i], 0, 0};
  }
}
} // namespace ELFF::Models
