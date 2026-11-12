#include "elff/models/rigidbody/RigidContactBody.hpp"
#include <cmath>

namespace ELFF::Models {
std::vector<Contact::ShapeRecord> RigidContactBody::contact_shapes() const {
  const auto* provider = dynamic_cast<const Contact::CollisionGeometryProvider*>(this);
  ELFF_VERIFY(provider, "Rigid contact participant requires collision geometry.");
  auto shape = provider->collision_shape();
  const auto& body = rigid_body();
  std::visit(
      [&](auto& geometry) {
        using Shape = std::decay_t<decltype(geometry)>;
        if constexpr (std::is_same_v<Shape, Contact::SphereGeometry>)
          geometry.center = body.world_point(geometry.center);
        else {
          geometry.axis.start = body.world_point(geometry.axis.start);
          geometry.axis.end = body.world_point(geometry.axis.end);
        }
      },
      shape);
  return {{0, shape, body.center_of_mass(), body.com_velocity(), body.angular_velocity_world()}};
}
void RigidContactBody::add_contact_force(size_t shape, real_t, const Contact::Vec3& point,
                                         const Contact::Vec3& force) {
  ELFF_VERIFY(shape == 0, "Invalid rigid shape ID.");
  const auto* provider = dynamic_cast<const Contact::CollisionGeometryProvider*>(this);
  const auto geometry = provider->collision_shape();
  const auto* sphere = std::get_if<Contact::SphereGeometry>(&geometry);
  const bool centered = sphere && sphere->center == Contact::Vec3{};
  Contact::Vec3 lever;
  for (size_t d = 0; d < 3; ++d) {
    lever[d] = point[d] - rigid_body().center_of_mass()[d];
    load_.force[d] += force[d];
  }
  if (!centered)
    for (size_t d = 0; d < 3; ++d) {
      const size_t e = (d + 1) % 3, f = (d + 2) % 3;
      load_.torque[d] += lever[e] * force[f] - lever[f] * force[e];
    }
  for (size_t d = 0; d < 3; ++d)
    ELFF_VERIFY(std::isfinite(load_.force[d]) && std::isfinite(load_.torque[d]),
                "Rigid contact load overflow.");
}
void RigidContactBody::predict_contact(real_t dt, const Contact::Vec3& gravity) {
  auto load = load_;
  for (size_t d = 0; d < 3; ++d)
    load.force[d] += rigid_body().mass() * gravity[d];
  rigid_body().kick(dt / 2, load);
  rigid_body().drift(dt);
}
void RigidContactBody::finish_contact(real_t dt, const Contact::Vec3& gravity) {
  auto load = load_;
  for (size_t d = 0; d < 3; ++d)
    load.force[d] += rigid_body().mass() * gravity[d];
  rigid_body().kick(dt / 2, load);
}
void RigidContactBody::advance_contact(real_t dt, const Contact::Vec3& gravity) {
  predict_contact(dt, gravity);
  finish_contact(dt, gravity);
}
void RigidContactBody::synchronize_contact_time(real_t time) {
  ELFF_VERIFY(std::isfinite(time) && time >= 0, "Invalid contact clock.");
  rigid_body().set_time_integration_state(time, contact_iteration());
}
void RigidContactBody::validate_contact_participation() const {
  ELFF_VERIFY(dynamic_cast<const Contact::CollisionGeometryProvider*>(this),
              "Missing collision geometry.");
  rigid_body().verify_axisymmetric_inertia();
}
} // namespace ELFF::Models
