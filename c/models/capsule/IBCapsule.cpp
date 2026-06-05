#include "elff/c/models/capsule/IBCapsule.h"

#include "elff/config/config.hpp"
#include "elff/general/error.hpp"
#include "elff/models/capsule/IBMCapsule.hpp"

#include <memory>

using namespace ELFF;
using namespace ELFF::Models;

namespace {

Vec3 to_vec3(vertex_t value)
{
  return { static_cast<double>(value.x),
           static_cast<double>(value.y),
           static_cast<double>(value.z) };
}

IBMCapsule* to_capsule(ib_capsule_t handle)
{
  return reinterpret_cast<IBMCapsule*>(handle);
}

void verify_common_inputs(int refinements, double elastic_modulus)
{
  ELFF_VERIFY(refinements >= 0,
              "ib_capsule_*_new(): refinements must be nonnegative.\n");
  ELFF_VERIFY(elastic_modulus >= 0.0,
              "ib_capsule_*_new(): elastic_modulus must be nonnegative.\n");
}

void reinitialize_structure(IBMCapsule* capsule)
{
  ELFF_VERIFY(capsule != nullptr, "ib_capsule: null capsule handle.\n");
  capsule->initializeStructuralState();
}

} // namespace

extern "C"
{

ib_capsule_t
ib_capsule_sphere_new(double radius,
                      vertex_t center,
                      int refinements,
                      double elastic_modulus)
{
  verify_common_inputs(refinements, elastic_modulus);
  ELFF_VERIFY(radius > 0.0,
              "ib_capsule_sphere_new(): radius must be positive.\n");

  CapsuleSphereSpec spec;
  spec.radius = radius;
  spec.center = to_vec3(center);
  spec.refinements = refinements;

  return new IBMCapsule(
      CapsuleMeshBuilder::sphere(spec),
      std::make_unique<CapsuleNeoHookeanLaw>(elastic_modulus));
}

ib_capsule_t
ib_capsule_ellipsoid_new(vertex_t radii,
                         vertex_t center,
                         int refinements,
                         double elastic_modulus)
{
  verify_common_inputs(refinements, elastic_modulus);
  CapsuleEllipsoidSpec spec;
  spec.radii = to_vec3(radii);
  spec.center = to_vec3(center);
  spec.refinements = refinements;

  return new IBMCapsule(
      CapsuleMeshBuilder::ellipsoid(spec),
      std::make_unique<CapsuleNeoHookeanLaw>(elastic_modulus));
}

ib_capsule_t
ib_capsule_biconcave_new(double radius,
                         vertex_t center,
                         int refinements,
                         double elastic_modulus)
{
  verify_common_inputs(refinements, elastic_modulus);
  ELFF_VERIFY(radius > 0.0,
              "ib_capsule_biconcave_new(): radius must be positive.\n");

  CapsuleBiconcaveSpec spec;
  spec.radius = radius;
  spec.center = to_vec3(center);
  spec.refinements = refinements;

  return new IBMCapsule(
      CapsuleMeshBuilder::biconcave(spec),
      std::make_unique<CapsuleNeoHookeanLaw>(elastic_modulus));
}

void
ib_capsule_set_neo_hookean_law(ib_capsule_t handle, double elastic_modulus)
{
  ELFF_VERIFY(elastic_modulus >= 0.0,
              "ib_capsule_set_neo_hookean_law(): elastic_modulus must be "
              "nonnegative.\n");
  IBMCapsule* capsule = to_capsule(handle);
  capsule->setLaw(std::make_unique<CapsuleNeoHookeanLaw>(elastic_modulus));
  reinitialize_structure(capsule);
}

void
ib_capsule_set_skalak_law(ib_capsule_t handle,
                          double elastic_modulus,
                          double area_dilatation_modulus)
{
  ELFF_VERIFY(elastic_modulus >= 0.0,
              "ib_capsule_set_skalak_law(): elastic_modulus must be "
              "nonnegative.\n");
  ELFF_VERIFY(area_dilatation_modulus >= 0.0,
              "ib_capsule_set_skalak_law(): area_dilatation_modulus must be "
              "nonnegative.\n");
  IBMCapsule* capsule = to_capsule(handle);
  capsule->setLaw(std::make_unique<CapsuleSkalakLaw>(
      elastic_modulus, area_dilatation_modulus));
  reinitialize_structure(capsule);
}

void
ib_capsule_set_linear_bending_law(ib_capsule_t handle, double bending_modulus)
{
  ELFF_VERIFY(bending_modulus >= 0.0,
              "ib_capsule_set_linear_bending_law(): bending_modulus must be "
              "nonnegative.\n");
  IBMCapsule* capsule = to_capsule(handle);
  capsule->setBendingLaw(
      std::make_unique<CapsuleLinearBendingLaw>(bending_modulus));
  capsule->initializeReferenceCurvature();
}

void
ib_capsule_set_helfrich_bending_law(ib_capsule_t handle,
                                    double bending_modulus)
{
  ELFF_VERIFY(bending_modulus >= 0.0,
              "ib_capsule_set_helfrich_bending_law(): bending_modulus must be "
              "nonnegative.\n");
  IBMCapsule* capsule = to_capsule(handle);
  capsule->setBendingLaw(
      std::make_unique<CapsuleHelfrichBendingLaw>(bending_modulus));
  capsule->initializeReferenceCurvature();
}

void
ib_capsule_set_constant_reference_curvature(ib_capsule_t handle,
                                            double reference_curvature)
{
  IBMCapsule* capsule = to_capsule(handle);
  capsule->setConstantReferenceCurvature(reference_curvature);
}

void
ib_capsule_clear_bending_law(ib_capsule_t handle)
{
  to_capsule(handle)->setBendingLaw(nullptr);
}

int
ib_capsule_get_triangle_count(ib_capsule_t handle)
{
  IBMCapsule* capsule = to_capsule(handle);
  ELFF_VERIFY(capsule != nullptr,
              "ib_capsule_get_triangle_count(): null capsule handle.\n");
  return capsule->mesh().numTriangles();
}

void
ib_capsule_get_triangle_node_ids(ib_capsule_t handle,
                                 int triangle_index,
                                 int node_ids[3])
{
  IBMCapsule* capsule = to_capsule(handle);
  ELFF_VERIFY(capsule != nullptr,
              "ib_capsule_get_triangle_node_ids(): null capsule handle.\n");
  ELFF_VERIFY(node_ids != nullptr,
              "ib_capsule_get_triangle_node_ids(): null node_ids buffer.\n");
  ELFF_VERIFY(triangle_index >= 0 &&
                  triangle_index < capsule->mesh().numTriangles(),
              "ib_capsule_get_triangle_node_ids(): invalid triangle index.\n");

  const auto& triangle =
      capsule->mesh().triangles[static_cast<size_t>(triangle_index)];
  node_ids[0] = triangle.nodes[0];
  node_ids[1] = triangle.nodes[1];
  node_ids[2] = triangle.nodes[2];
}

ib_velocity_coupled_t
ib_capsule_as_velocity_coupled(ib_capsule_t handle)
{
  return static_cast<IBVelocityCoupled*>(to_capsule(handle));
}

void
ib_capsule_destroy(ib_capsule_t handle)
{
  delete to_capsule(handle);
}

} // extern "C"
