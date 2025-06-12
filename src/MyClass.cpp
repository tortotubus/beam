// src/MyClass.cpp
#include "beam/MyClass.hpp"

namespace beam {

MyClass::MyClass(int init_value)
  : _base(init_value)
{}

MyClass::~MyClass() = default;

int MyClass::compute(int x) const {
  // just an example: add base and square
  return _base + x*x;
}

} // namespace beam
