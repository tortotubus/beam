#include "cursor.hpp"
#include <algorithm>

namespace editor {

  void Cursor::clamp(const Projection &projection) {
  selected = projection.locations.empty() ? 0 : std::min(selected, projection.locations.size() - 1);
}

void Cursor::move(Motion motion, const Projection &projection) {
  clamp(projection);
  bool upward = motion == Motion::Up;
  bool downward = motion == Motion::Down;

  if (!upward && !downward)
    resetColumn();

  if (projection.locations.empty())
    return;

  const auto &locations = projection.locations;
  const auto &location = locations[selected];
  
  if (upward || downward) {
    if (!preferredColumn)
      preferredColumn = location.column;
    std::size_t bestLineDistance = std::string::npos;
    std::size_t bestColumnDistance = std::string::npos;
    for (std::size_t i = 0; i < locations.size(); ++i) {
      const auto &candidate = locations[i];
      if (upward ? candidate.line >= location.line
                 : candidate.line <= location.line)
        continue;
      auto lineDistance = upward ? location.line - candidate.line
                                 : candidate.line - location.line;
      auto columnDistance = std::max(candidate.column, *preferredColumn) -
                            std::min(candidate.column, *preferredColumn);
      if (lineDistance < bestLineDistance ||
          (lineDistance == bestLineDistance &&
           columnDistance < bestColumnDistance)) {
        selected = i;
        bestLineDistance = lineDistance;
        bestColumnDistance = columnDistance;
      }
    }
  } else if (motion == Motion::Left || motion == Motion::Right) {
    std::size_t bestDistance = std::string::npos;
    bool left = motion == Motion::Left;
    for (std::size_t i = 0; i < locations.size(); ++i) {
      const auto &candidate = locations[i];
      if (candidate.line != location.line)
        continue;
      if (left ? candidate.column >= location.column
               : candidate.column <= location.column)
        continue;
      auto distance = left ? location.column - candidate.column
                           : candidate.column - location.column;
      if (distance < bestDistance) {
        selected = i;
        bestDistance = distance;
      }
    }
  } else if (motion == Motion::Previous && selected)
    --selected;
  else if (motion == Motion::Next && selected + 1 < locations.size())
    ++selected;
  else if (motion == Motion::First)
    selected = 0;
  else if (motion == Motion::Last)
    selected = locations.size() - 1;
}
} // namespace editor
