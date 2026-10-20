#include "elff/c/models/rigidbody/IBRigidBody.h"

#include "elff/config/config.hpp"
#include "elff/general/error.hpp"
#include "elff/models/rigidbody/IBPinnedRigidBodyCircle.hpp"
#include "elff/models/rigidbody/IBPinnedRigidBodyFibre.hpp"
#include "elff/models/rigidbody/IBPinnedRigidBodySphere.hpp"
#include "elff/models/rigidbody/IBRigidBodyCircle.hpp"
#include "elff/models/rigidbody/IBRigidBodyFibre.hpp"
#include "elff/models/rigidbody/IBRigidBodySphere.hpp"
#include "elff/models/ibm/IBCoupling.hpp"

#include <memory>

using namespace ELFF;
using namespace ELFF::Models;

namespace {
RigidBodyModel::Vec3
to_vec3(vertex_t v)
{
  return { static_cast<real_t>(v.x),
           static_cast<real_t>(v.y),
           static_cast<real_t>(v.z) };
}

IBRigidBodySphere::Quat
to_quat(quaternion_t q)
{
  return { static_cast<real_t>(q.w),
           static_cast<real_t>(q.x),
           static_cast<real_t>(q.y),
           static_cast<real_t>(q.z) };
}
} // namespace

extern "C"
{

ib_rigid_body_circle_t
ib_rigid_body_circle_new(ib_coupling_t runtime,
                         double radius,
                         int point_count,
                         double density,
                         vertex_t center,
                         vertex_t velocity,
                         double angle,
                         double angular_velocity_z)
{
  if (!runtime)
    return nullptr;
  try {
  ELFF_VERIFY(point_count > 0,
              "ib_rigid_body_circle_new(): point_count must be positive.\n");

  auto body =
    std::make_unique<IBRigidBodyCircle>(static_cast<real_t>(radius),
                          static_cast<size_t>(point_count),
                          static_cast<real_t>(density),
                          to_vec3(center),
                          to_vec3(velocity),
                          static_cast<real_t>(angle),
                          static_cast<real_t>(angular_velocity_z));

  return reinterpret_cast<IBCoupling*>(runtime)->adopt_owned_model(
    ModelType::RigidBodyCircle, std::move(body));
  } catch (...) {
    return nullptr;
  }
}

ib_rigid_body_fibre_t
ib_rigid_body_fibre_new(ib_coupling_t runtime,
                        double length,
                        double diameter,
                        int nodes,
                        double linear_density,
                        vertex_t center,
                        vertex_t velocity,
                        quaternion_t q_bw,
                        vertex_t angular_velocity_body)
{
  if (!runtime)
    return nullptr;
  try {
  ELFF_VERIFY(nodes >= 2,
              "ib_rigid_body_fibre_new(): nodes must be at least 2.\n");

  auto body =
    std::make_unique<IBRigidBodyFibre>(static_cast<real_t>(length),
                         static_cast<real_t>(diameter),
                         static_cast<size_t>(nodes),
                         static_cast<real_t>(linear_density),
                         to_vec3(center),
                         to_vec3(velocity),
                         to_quat(q_bw),
                         to_vec3(angular_velocity_body));

  return reinterpret_cast<IBCoupling*>(runtime)->adopt_owned_model(
    ModelType::RigidBodyFibre, std::move(body));
  } catch (...) {
    return nullptr;
  }
}

ib_pinned_rigid_body_circle_t
ib_pinned_rigid_body_circle_new(ib_coupling_t runtime,
                                double radius,
                                int point_count,
                                double density,
                                vertex_t center,
                                double angle)
{
  if (!runtime)
    return nullptr;
  try {
  ELFF_VERIFY(point_count > 0,
              "ib_pinned_rigid_body_circle_new(): point_count must be "
              "positive.\n");

  auto body =
    std::make_unique<IBPinnedRigidBodyCircle>(static_cast<real_t>(radius),
                                static_cast<size_t>(point_count),
                                static_cast<real_t>(density),
                                to_vec3(center),
                                static_cast<real_t>(angle));

  return reinterpret_cast<IBCoupling*>(runtime)->adopt_owned_model(
    ModelType::PinnedRigidBodyCircle, std::move(body));
  } catch (...) {
    return nullptr;
  }
}

ib_pinned_rigid_body_fibre_t
ib_pinned_rigid_body_fibre_new(ib_coupling_t runtime,
                               double length,
                               double diameter,
                               int nodes,
                               double linear_density,
                               vertex_t center,
                               quaternion_t q_bw)
{
  if (!runtime)
    return nullptr;
  try {
  ELFF_VERIFY(nodes >= 2,
              "ib_pinned_rigid_body_fibre_new(): nodes must be at least 2.\n");

  auto body =
    std::make_unique<IBPinnedRigidBodyFibre>(static_cast<real_t>(length),
                               static_cast<real_t>(diameter),
                               static_cast<size_t>(nodes),
                               static_cast<real_t>(linear_density),
                               to_vec3(center),
                               to_quat(q_bw));

  return reinterpret_cast<IBCoupling*>(runtime)->adopt_owned_model(
    ModelType::PinnedRigidBodyFibre, std::move(body));
  } catch (...) {
    return nullptr;
  }
}

ib_rigid_body_sphere_t
ib_rigid_body_sphere_new(ib_coupling_t runtime,
                         double radius,
                         int point_count,
                         double density,
                         vertex_t center,
                         vertex_t velocity,
                         quaternion_t q_bw,
                         vertex_t angular_velocity_body)
{
  if (!runtime)
    return nullptr;
  try {
  ELFF_VERIFY(point_count > 0,
              "ib_rigid_body_sphere_new(): point_count must be positive.\n");

  auto body =
    std::make_unique<IBRigidBodySphere>(static_cast<real_t>(radius),
                          static_cast<size_t>(point_count),
                          static_cast<real_t>(density),
                          to_vec3(center),
                          to_vec3(velocity),
                          to_quat(q_bw),
                          to_vec3(angular_velocity_body));

  return reinterpret_cast<IBCoupling*>(runtime)->adopt_owned_model(
    ModelType::RigidBodySphere, std::move(body));
  } catch (...) {
    return nullptr;
  }
}

ib_pinned_rigid_body_sphere_t
ib_pinned_rigid_body_sphere_new(ib_coupling_t runtime,
                                double radius,
                                int point_count,
                                double density,
                                vertex_t center,
                                quaternion_t q_bw)
{
  if (!runtime)
    return nullptr;
  try {
  ELFF_VERIFY(point_count > 0,
              "ib_pinned_rigid_body_sphere_new(): point_count must be "
              "positive.\n");

  auto body =
    std::make_unique<IBPinnedRigidBodySphere>(static_cast<real_t>(radius),
                                static_cast<size_t>(point_count),
                                static_cast<real_t>(density),
                                to_vec3(center),
                                to_quat(q_bw));

  return reinterpret_cast<IBCoupling*>(runtime)->adopt_owned_model(
    ModelType::PinnedRigidBodySphere, std::move(body));
  } catch (...) {
    return nullptr;
  }
}

} // extern "C"
