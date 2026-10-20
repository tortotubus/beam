#pragma once

#include "elff/models/Model.hpp"

namespace ELFF::Models {
/** @brief Marker interface for runtime-owned immersed-boundary models. */
class IBModel : public Model {
public:
  ~IBModel() override = default;
};
} // namespace ELFF::Models
