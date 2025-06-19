#pragma once

#include <iostream>
#include <limits>

namespace beam {
inline void skip_comment_lines(std::istream &is, const char comment_char) 
{
  while (1) {
    is >> std::ws;
    if (is.peek() != comment_char) {
      break;
    }
    is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  }
}
} // namespace beam