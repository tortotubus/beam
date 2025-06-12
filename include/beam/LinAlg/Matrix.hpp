// Matrix.hpp
#pragma once

#include <cstddef>
#include <iostream>
#include <stdexcept>
#include <utility>

namespace beam {

enum MatrixOrder {
  RowMajorOrder = 0,
  ColMajorOrder = 1,
};

template <class T> class Matrix {
public:
  Matrix(std::size_t rows, std::size_t cols, MatrixOrder order);

  // Initializer list constructor
  Matrix(std::initializer_list<std::initializer_list<T>> init_list);

  // Copy Constructor (Deep Copy)
  Matrix(const Matrix<T> &other);

  // Copy Assignment Operator (Deep Copy)
  Matrix &operator=(const Matrix<T> &other);

  // Move Constructor
  Matrix(Matrix<T> &&other) noexcept;

  // Move Assignment Operator
  Matrix &operator=(Matrix<T> &&other) noexcept;

  ~Matrix();

  // element access
  T &operator()(std::size_t row, std::size_t col);
  const T &operator()(std::size_t row, std::size_t col) const;

  // dimensions
  std::size_t size() const;
  std::size_t rows() const;
  std::size_t cols() const;

  // pretty‐print
  friend std::ostream &operator<<(std::ostream &os, const Matrix &m) {
    for (std::size_t r = 0; r < m.nrows; ++r) {
      for (std::size_t c = 0; c < m.ncols; ++c) {
        os << m(r, c) << "\t";
      }
      os << "\n";
    }
    return os;
  }

private:
  std::size_t nrows, ncols, n;
  MatrixOrder order;
  T *data;

  std::size_t index(std::size_t row, std::size_t col) const;
  bool in_bounds(std::size_t row, std::size_t col) const;
};

} // namespace beam
