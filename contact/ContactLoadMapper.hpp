#pragma once

namespace ELFF::Contact {

/**
 * @brief Maps contact forces and moments to model loads.
 *
 * Implementations should transfer contact forces and moments to the owning
 * model's nodal or generalized loads while preserving virtual work. Contact
 * loads are integrated mechanical forces and must not receive the quadrature
 * weights used to convert immersed-boundary force densities.
 */
class ContactLoadMapper {
public:
  /**
   * @brief Destroy the interface polymorphically.
   */
  virtual ~ContactLoadMapper() = 0;
};

inline ContactLoadMapper::~ContactLoadMapper() = default;

} // namespace ELFF::Contact
