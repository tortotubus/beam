#pragma once
#include "elff/config/config.hpp"
#include <array>

namespace ELFF::Contact {
using Vec3 = std::array<real_t, 3>;
/** @brief Integrated world force and torque, about the body's center of mass. */
struct Wrench {
  Vec3 force{};
  Vec3 torque{};
};
} // namespace ELFF::Contact
