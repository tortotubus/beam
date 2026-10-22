#pragma once
#include "projection.hpp"
#include <optional>

namespace editor {
/** @brief Spatial navigation over a projection, independent of terminal keys.
 */
class Cursor {
public:
  enum class Motion { Left, Right, Up, Down, Next, Previous, First, Last };
  /**
   * @brief Return the selected index in the current projection location list.
   */
  std::size_t index() const { return selected; }
  /**
   * @brief Keep the selected index in range after the projection changes.
   */
  void clamp(const Projection &projection);
  /**
   * @brief Move between locations, retaining the target column for vertical
   * motion.
   */
  void move(Motion motion, const Projection &projection);
  /**
   * @brief Forget the preferred column before a new vertical navigation
   * sequence.
   */
  void resetColumn() { preferredColumn.reset(); }

private:
  std::size_t selected = 0;
  std::optional<std::size_t> preferredColumn;
};
} // namespace editor
