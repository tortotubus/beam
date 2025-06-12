#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <type_traits>
#include <utility>

#include <iomanip> // for scientific, setprecision, setw, copyfmt
#include <iostream>

namespace beam {

template <class T> class Vector {
public:
  Vector(size_t n);

  // Constructor that accepts an initializer list
  Vector(std::initializer_list<T> init_list);

  // Copy Constructor (Deep Copy)
  Vector(const Vector &other);

  // Copy Assignment Operator (Deep Copy)
  Vector &operator=(const Vector &other);

  // Move Constructor
  Vector(Vector &&other) noexcept;

  // Move Assignment Operator
  Vector &operator=(Vector &&other) noexcept;

  // Destructor
  ~Vector();

  /*
      Access operators
  */
  T &operator()(size_t i);
  const T &operator()(size_t i) const;
  size_t size() const;

  // Add friend declaration for operator<<
  // inline friend definition
  friend std::ostream &operator<<(std::ostream &os, Vector const &v) {
    // 1) make a bare-bones ios object and copy the current state into it
    std::ios oldFmt(nullptr);
    oldFmt.copyfmt(os);

    // 2) set your new formatting
    os << std::scientific << std::setprecision(6);

    // 3) output with a fixed width
    os << "[";
    for (size_t i = 0; i < v.n; ++i) {
      os << std::setw(12) << v.data[i];
      if (i + 1 < v.n)
        os << ", ";
    }
    os << "]";

    // 4) restore the old formatting
    os.copyfmt(oldFmt);
    return os;
  }

private:
  size_t n;
  T *data;

  inline bool in_bounds(size_t i) const;
};

} // namespace beam