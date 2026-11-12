#include "elff/models/articulated/DryBendingRodImplicit.hpp"
#include <algorithm>
#include <cmath>
#include <limits>

namespace ELFF::Models {
DryBendingRodImplicit::DryBendingRodImplicit(const BendingRodImplicitConfig& config, real_t radius,
                                             bool body_contacts)
    : BendingRodImplicit(config), radius_(radius), body_contacts_(body_contacts),
      loads_(get_mesh().get_nodes()) {
  validate_contact_participation();
}
void DryBendingRodImplicit::validate_contact_participation() const {
  ELFF_VERIFY(std::isfinite(radius_) && radius_ > 0,
              "Collision radius must be finite and positive.");
}
std::vector<Contact::ShapeRecord> DryBendingRodImplicit::contact_shapes() const {
  const auto& mesh = get_mesh();
  std::vector<Contact::ShapeRecord> records;
  for (size_t i = 0; i + 1 < mesh.get_nodes(); ++i) {
    const auto& a = mesh.get_centerline(i);
    const auto& b = mesh.get_centerline(i + 1);
    Contact::Vec3 middle;
    for (size_t d = 0; d < 3; ++d)
      middle[d] = a[d] + (b[d] - a[d]) / 2;
    Contact::ShapeRecord shape{0, Contact::CapsuleGeometry{{a, b}, radius_}, middle, {}, {}};
    shape.local_shape = i;
    shape.endpoint_velocities = std::array<Contact::Vec3, 2>{mesh.get_centerline_velocity(i),
                                                             mesh.get_centerline_velocity(i + 1)};
    shape.endpoint_features = std::array<size_t, 2>{i, i + 1};
    shape.body_contacts = body_contacts_;
    records.push_back(shape);
  }
  return records;
}
void DryBendingRodImplicit::clear_contact_loads() {
  std::fill(loads_.begin(), loads_.end(), Vector{});
}
void DryBendingRodImplicit::add_contact_force(size_t shape, real_t coordinate, const Contact::Vec3&,
                                              const Contact::Vec3& force) {
  ELFF_VERIFY(shape + 1 > shape && shape + 1 < loads_.size() && std::isfinite(coordinate) &&
                  coordinate >= 0 && coordinate <= 1,
              "Invalid segment contact coordinate.");
  for (size_t d = 0; d < 3; ++d) {
    ELFF_VERIFY(std::isfinite(force[d]), "Contact force must be finite.");
    loads_[shape][d] += (1 - coordinate) * force[d];
    loads_[shape + 1][d] += coordinate * force[d];
    ELFF_VERIFY(std::isfinite(loads_[shape][d]) && std::isfinite(loads_[shape + 1][d]),
                "Nodal load overflow.");
  }
}
Contact::Wrench DryBendingRodImplicit::contact_wrench() const {
  Contact::Wrench result;
  const auto& origin = get_mesh().get_centerline(0);
  for (size_t i = 0; i < loads_.size(); ++i) {
    Contact::Vec3 lever;
    for (size_t d = 0; d < 3; ++d) {
      lever[d] = get_mesh().get_centerline(i)[d] - origin[d];
      result.force[d] += loads_[i][d];
    }
    for (size_t d = 0; d < 3; ++d) {
      const size_t e = (d + 1) % 3, f = (d + 2) % 3;
      result.torque[d] += lever[e] * loads_[i][f] - lever[f] * loads_[i][e];
    }
  }
  return result;
}
void DryBendingRodImplicit::validate_step(real_t dt) const {
  ELFF_VERIFY(
      std::isfinite(dt) && dt >= 0 && std::isfinite(time_ + dt) &&
          (dt == 0 || (time_ + dt > time_ && iteration_ < std::numeric_limits<size_t>::max())),
      "Invalid rod timestep.");
}
void DryBendingRodImplicit::solve(real_t dt, const std::vector<Vector>& density) {
  validate_step(dt);
  BendingRodImplicit::solve(dt, density);
  if (dt > 0) {
    time_ += dt;
    ++iteration_;
  }
}
void DryBendingRodImplicit::solve_integrated(real_t dt, const std::vector<Vector>& forces) {
  validate_step(dt);
  BendingRodImplicit::solve_integrated(dt, forces);
  if (dt > 0) {
    time_ += dt;
    ++iteration_;
  }
}
void DryBendingRodImplicit::apply_initial_condition(const BendingRodMesh& mesh) {
  BendingRodImplicit::apply_initial_condition(mesh);
  time_ = 0;
  iteration_ = 0;
  clear_contact_loads();
}
void DryBendingRodImplicit::advance_contact(real_t dt, const Contact::Vec3& gravity) {
  auto forces = loads_;
  const auto masses = nodal_masses();
  for (size_t i = 0; i < forces.size(); ++i)
    for (size_t d = 0; d < 3; ++d)
      forces[i][d] += masses[i] * gravity[d];
  solve_integrated(dt, forces);
}
void DryBendingRodImplicit::synchronize_contact_time(real_t time) {
  ELFF_VERIFY(std::isfinite(time) && time >= 0, "Invalid rod clock.");
  time_ = time;
}
void DryBendingRodImplicit::pack_state(ModelState& state) const {
  ELFF_VERIFY(iteration_ <= static_cast<size_t>(std::numeric_limits<int64_t>::max()),
              "Rod clock overflow.");
  BendingRodImplicit::pack_state(state);
  state.ints.push_back(2);
  state.ints.push_back(static_cast<int64_t>(iteration_));
  state.ints.push_back(body_contacts_ ? 1 : 0);
  state.reals.push_back(radius_);
  state.reals.push_back(time_);
}
void DryBendingRodImplicit::unpack_state(const ModelState& state) {
  ELFF_VERIFY(state.ints.size() == 8 && state.ints[5] == 2 && state.ints[6] >= 0 &&
                  static_cast<uint64_t>(state.ints[6]) <= std::numeric_limits<size_t>::max() &&
                  state.ints[7] == (body_contacts_ ? 1 : 0) && state.reals.size() >= 2,
              "Invalid contact rod checkpoint.");
  const auto count = state.reals.size();
  ELFF_VERIFY(state.reals[count - 2] == radius_ && std::isfinite(state.reals.back()) &&
                  state.reals.back() >= 0,
              "Contact radius or clock mismatch.");
  ModelState structural = state;
  structural.ints.resize(5);
  structural.reals.resize(count - 2);
  BendingRodImplicit::unpack_state(structural);
  iteration_ = static_cast<size_t>(state.ints[6]);
  time_ = state.reals.back();
  clear_contact_loads();
}
void DryBendingRodImplicit::append_vtk_polydata(std::vector<IO::CXX::vtkPolyData>& datasets) const {
  IO::CXX::vtkPolyData data;
  for (const auto& p : get_mesh().get_centerline())
    data.add_point(static_cast<float>(p[0]), static_cast<float>(p[1]), static_cast<float>(p[2]));
  for (size_t i = 1; i < get_mesh().get_nodes(); ++i)
    data.add_line(static_cast<int64_t>(i - 1), static_cast<int64_t>(i));
  datasets.push_back(std::move(data));
}
} // namespace ELFF::Models
