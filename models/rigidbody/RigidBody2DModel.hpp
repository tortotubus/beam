#pragma once

#include "elff/models/rigidbody/RigidBody2DMesh.hpp"

namespace ELFF {
namespace Models {

/**
 * @brief Placeholder 2D rigid-body model API.
 *
 * TODO: Define 2D rigid-body state, Newton-Euler planar dynamics, and
 * time integration methods.
 */
class RigidBody2DModel
{
public:
  RigidBody2DModel() = default;
  virtual ~RigidBody2DModel() = default;

protected:
  RigidBody2DMesh mesh_;
};

} // namespace Models
} // namespace ELFF
