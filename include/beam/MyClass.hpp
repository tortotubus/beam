// include/mylib/MyClass.hpp
#pragma once

namespace beam {

class MyClass {
public:
  explicit MyClass(int init_value);
  ~MyClass();

  // do some work and return a result
  int compute(int x) const;

private:
  int _base;
};

} // namespace beam
