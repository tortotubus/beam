#include "dummy.hpp"

namespace beam {

Dummy::Dummy(int init_value) : _base(init_value) {}

Dummy::~Dummy() = default;

int Dummy::compute(int x) const {
  // just an example: add base and square
  return _base + x*x;
}

} // namespace beam
