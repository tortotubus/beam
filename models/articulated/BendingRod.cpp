#include "elff/models/articulated/BendingRod.hpp"

namespace ELFF::Models {
void BendingRod::solve(real_t dt) { solve(dt, Vector{0, 0, 0}); }
void BendingRod::solve(real_t dt, Vector load) {
  solve(dt, std::vector<Vector>(mesh_.get_nodes(), load));
}
} // namespace ELFF::Models
