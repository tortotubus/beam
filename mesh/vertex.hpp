#pragma once

#include "../config/config.hpp"

namespace beam {
namespace fem {
class Vertex
{
protected:
  real_t coord[3];

public:
  Vertex() = default;
  // Vertex(real_t *xx, int dim);
  Vertex(real_t x, real_t y)
  {
    coord[0] = x;
    coord[1] = y;
    coord[2] = 0.;
  }
  Vertex(real_t x, real_t y, real_t z)
  {
    coord[0] = x;
    coord[1] = y;
    coord[2] = z;
  }
  /// Returns pointer to the coordinates of the vertex.
  inline real_t* operator()() const { return (real_t*)coord; }
  /// Returns the i'th coordinate of the vertex.
  inline real_t& operator()(int i) { return coord[i]; }
  /// Returns the i'th coordinate of the vertex.
  inline const real_t& operator()(int i) const { return coord[i]; }
};
}
} // namespace beam
