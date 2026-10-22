#pragma once
#include "elff/contact/ContactGeometryProvider.hpp"

#include "elff/models/articulated/BendingRodMesh.hpp"

namespace ELFF::Models {
/** Common dynamic rod interface. Numerical coordinates and mechanics belong
 * to concrete implementations. Loads are force densities per unit length.
 * Edit a mesh and call apply_initial_condition to change the solver state;
 * editing the published mesh alone does not advance or reset the solver.
 */
class BendingRod : public Contact::ContactGeometryProvider {
public:
  /** @brief Borrow the current nodal polyline; reacquire after mesh updates. */
  Contact::PolylineView contact_geometry() const override {
    const auto& points = get_mesh().get_centerline();
    return {points.data(), points.size()};
  }
  using Vector = BendingRodMesh::Vector;
  enum class Boundary { Pinned, Free };
  virtual ~BendingRod() = default;
  virtual void solve(real_t dt);
  virtual void solve(real_t dt, Vector uniform_load);
  virtual void solve(real_t dt, const std::vector<Vector>& nonuniform_load) = 0;
  virtual void apply_initial_condition(const BendingRodMesh& mesh) = 0;
  virtual BendingRodMesh& get_mesh() { return mesh_; }
  virtual const BendingRodMesh& get_mesh() const { return mesh_; }

protected:
  BendingRod(size_t nodes, real_t length) : mesh_(nodes, length) {}
  BendingRodMesh mesh_;
};
} // namespace ELFF::Models
