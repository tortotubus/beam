#include "elff/contact/ContactSystem.hpp"
#include "elff/general/error.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <tuple>

namespace ELFF::Contact {

ContactSystem::ContactSystem(const ContactGrid& grid, ContactParameters parameters)
    : parameters_(parameters), broad_phase_(grid.origin, grid.cell_width, grid.dimensions) {
  ELFF_VERIFY(std::isfinite(parameters.normal_stiffness) && parameters.normal_stiffness > 0 &&
                  std::isfinite(parameters.normal_damping) && parameters.normal_damping >= 0,
              "Contact stiffness must be positive and damping nonnegative; both finite.");
}

void ContactSystem::evaluate(const std::vector<ShapeRecord>& bodies,
                             const std::vector<PlaneGeometry>& planes) {
  contacts_.clear();
  shape_wrenches_.assign(bodies.size(), {});
  plane_reactions_.assign(planes.size(), {});
  broad_phase_.clear();
  std::set<std::pair<ShapeId, size_t>> ids;
  for (size_t i = 0; i < bodies.size(); ++i) {
    ELFF_VERIFY(ids.emplace(bodies[i].id, bodies[i].local_shape).second,
                "Owner/local shape IDs must be unique.");
    for (size_t d = 0; d < 3; ++d)
      ELFF_VERIFY(std::isfinite(bodies[i].center[d]) && std::isfinite(bodies[i].velocity[d]) &&
                      std::isfinite(bodies[i].angular_velocity[d]),
                  "Body kinematics must be finite.");
    if (bodies[i].endpoint_velocities) {
      ELFF_VERIFY(std::holds_alternative<CapsuleGeometry>(bodies[i].geometry),
                  "Endpoint kinematics require capsule geometry.");
      for (const auto& velocity : *bodies[i].endpoint_velocities)
        for (auto value : velocity)
          ELFF_VERIFY(std::isfinite(value), "Invalid endpoint velocity.");
    }
    broad_phase_.insert(i, collision_bounds(bodies[i].geometry));
  }
  auto velocity_at = [&](size_t i, const Vec3& point, real_t coordinate) {
    if (bodies[i].endpoint_velocities) {
      Vec3 velocity;
      for (size_t d = 0; d < 3; ++d)
        velocity[d] = (1 - coordinate) * (*bodies[i].endpoint_velocities)[0][d] +
                      coordinate * (*bodies[i].endpoint_velocities)[1][d];
      return velocity;
    }
    Vec3 lever, velocity = bodies[i].velocity;
    for (size_t d = 0; d < 3; ++d)
      lever[d] = point[d] - bodies[i].center[d];
    for (size_t d = 0; d < 3; ++d) {
      const size_t e = (d + 1) % 3, f = (d + 2) % 3;
      velocity[d] +=
          bodies[i].angular_velocity[e] * lever[f] - bodies[i].angular_velocity[f] * lever[e];
    }
    return velocity;
  };
  auto add_load = [](Wrench& load, const Vec3& center, const Vec3& point, const Vec3& force,
                     real_t sign, bool zero_torque) {
    Vec3 lever;
    for (size_t d = 0; d < 3; ++d) {
      lever[d] = point[d] - center[d];
      load.force[d] += sign * force[d];
    }
    if (!zero_torque)
      for (size_t d = 0; d < 3; ++d) {
        const size_t e = (d + 1) % 3, f = (d + 2) % 3;
        load.torque[d] += sign * (lever[e] * force[f] - lever[f] * force[e]);
      }
  };
  auto centered_sphere = [&](size_t i) {
    const auto* sphere = std::get_if<SphereGeometry>(&bodies[i].geometry);
    return sphere && sphere->center == bodies[i].center;
  };
  auto assemble = [&](size_t i, size_t j, bool plane, const ContactPoint& point) {
    if (point.gap >= 0)
      return;
    const auto va = velocity_at(i, point.point, point.parameter_a);
    const auto vb = plane ? Vec3{} : velocity_at(j, point.point, point.parameter_b);
    real_t speed = 0;
    for (size_t d = 0; d < 3; ++d)
      speed += (va[d] - vb[d]) * point.normal[d];
    const real_t trial =
        -parameters_.normal_stiffness * point.gap - parameters_.normal_damping * speed;
    ELFF_VERIFY(std::isfinite(trial), "Contact force overflow.");
    const real_t magnitude = std::max(real_t{0}, trial);
    Vec3 force;
    for (size_t d = 0; d < 3; ++d)
      force[d] = magnitude * point.normal[d];
    add_load(shape_wrenches_[i], bodies[i].center, point.point, force, 1, centered_sphere(i));
    if (plane)
      add_load(plane_reactions_[j], planes[j].point, point.point, force, -1, false);
    else
      add_load(shape_wrenches_[j], bodies[j].center, point.point, force, -1, centered_sphere(j));
    contacts_.push_back({bodies[i].id, plane ? static_cast<ShapeId>(j) : bodies[j].id, plane, point,
                         force, bodies[i].local_shape, plane ? 0 : bodies[j].local_shape});
  };
  // Canonical owner/shape ordering also fixes the sign convention for identical
  // static overlapping capsules. Group contacts by material endpoint features.
  using Feature = std::pair<bool, size_t>;
  using PairKey = std::tuple<ShapeId, Feature, ShapeId, Feature>;
  using Candidate = std::tuple<size_t, size_t, ContactPoint>;
  std::map<PairKey, Candidate> pairs;
  auto feature = [](const ShapeRecord& shape, real_t coordinate) -> Feature {
    const real_t tolerance = 64 * std::numeric_limits<real_t>::epsilon();
    if (shape.endpoint_features && coordinate <= tolerance)
      return {true, (*shape.endpoint_features)[0]};
    if (shape.endpoint_features && coordinate >= 1 - tolerance)
      return {true, (*shape.endpoint_features)[1]};
    return {false, shape.local_shape};
  };
  const auto candidates = broad_phase_.candidate_pairs_if([&](ShapeId i, ShapeId j) {
    return bodies[i].id != bodies[j].id && bodies[i].body_contacts && bodies[j].body_contacts;
  });
  for (const auto& candidate : candidates) {
    size_t i = candidate.first, j = candidate.second;
    if (std::tie(bodies[j].id, bodies[j].local_shape) <
        std::tie(bodies[i].id, bodies[i].local_shape))
      std::swap(i, j);
    const auto& a = bodies[i];
    const auto& b = bodies[j];
    auto point = collision_pair(a.geometry, b.geometry);
    if (point.gap >= 0)
      continue;
    if (point.normal_fallback) {
      const auto va = velocity_at(i, point.point, point.parameter_a);
      const auto vb = velocity_at(j, point.point, point.parameter_b);
      Vec3 hint;
      for (size_t d = 0; d < 3; ++d)
        hint[d] = vb[d] - va[d];
      // Oppose incoming relative motion when geometry alone cannot choose a side.
      point = collision_pair(a.geometry, b.geometry, hint);
    }
    const PairKey key{a.id, feature(a, point.parameter_a), b.id, feature(b, point.parameter_b)};
    auto found = pairs.find(key);
    bool replace = found == pairs.end();
    if (!replace) {
      const auto& previous = std::get<2>(found->second);
      const auto pi = std::get<0>(found->second), pj = std::get<1>(found->second);
      replace = point.gap < previous.gap ||
                (point.gap == previous.gap &&
                 std::tie(a.local_shape, b.local_shape) <
                     std::tie(bodies[pi].local_shape, bodies[pj].local_shape));
    }
    if (replace)
      pairs.insert_or_assign(key, Candidate{i, j, point});
  }
  for (const auto& pair : pairs)
    assemble(std::get<0>(pair.second), std::get<1>(pair.second), false, std::get<2>(pair.second));
  for (size_t j = 0; j < planes.size(); ++j) {
    const auto plane = make_plane(planes[j].point, planes[j].normal);
    std::map<std::pair<ShapeId, size_t>, std::pair<size_t, ContactPoint>> joints;
    for (size_t i = 0; i < bodies.size(); ++i) {
      const auto point = collision_plane(bodies[i].geometry, plane);
      if (point.gap >= 0)
        continue;
      const auto& features = bodies[i].endpoint_features;
      if (features && (point.parameter_a == 0 || point.parameter_a == 1)) {
        const auto key = std::make_pair(bodies[i].id, (*features)[point.parameter_a == 0 ? 0 : 1]);
        auto found = joints.find(key);
        if (found == joints.end() || point.gap < found->second.second.gap ||
            (point.gap == found->second.second.gap &&
             bodies[i].local_shape < bodies[found->second.first].local_shape))
          joints.insert_or_assign(key, std::make_pair(i, point));
      } else
        assemble(i, j, true, point);
    }
    for (const auto& entry : joints)
      assemble(entry.second.first, j, true, entry.second.second);
  }
  for (const auto* loads : {&shape_wrenches_, &plane_reactions_})
    for (const auto& load : *loads)
      for (size_t d = 0; d < 3; ++d)
        ELFF_VERIFY(std::isfinite(load.force[d]) && std::isfinite(load.torque[d]),
                    "Contact load overflow.");
}

} // namespace ELFF::Contact
