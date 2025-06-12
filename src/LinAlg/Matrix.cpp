// Matrix.cpp

#include <beam/LinAlg/Matrix.hpp>
#include <algorithm>
#include <stdexcept>
#include <initializer_list>
#include <sstream>

namespace beam {

template <class T>
Matrix<T>::Matrix(std::size_t rows, std::size_t cols, MatrixOrder order)
    : nrows(rows), ncols(cols), n(rows * cols), order(order), data(new T[n]()) {}

template <class T>
Matrix<T>::Matrix(std::initializer_list<std::initializer_list<T>> init_list)
    : nrows(init_list.size()),
      ncols(init_list.begin()->size()),
      n(init_list.size() * init_list.begin()->size()),
      data(new T[n]) {
  // Validate that all rows have the same number of columns
  for (const auto &row : init_list) {
    if (row.size() != ncols) {
      throw std::invalid_argument("Matrix initializer list must have rows of the same size.");
    }
  }
  // Flatten the list
  size_t idx = 0;
  for (const auto &row : init_list) {
    for (const auto &elem : row) {
      data[idx++] = elem;
    }
  }
  order = RowMajorOrder;
}

template <class T>
Matrix<T>::Matrix(const Matrix &other)
    : nrows(other.nrows), ncols(other.ncols), n(other.n), order(other.order) {
  data = new T[n];
  std::copy(other.data, other.data + n, data);
}

template <class T>
Matrix<T> &Matrix<T>::operator=(const Matrix &other) {
  if (this == &other)
    return *this;
  delete[] data;
  nrows = other.nrows;
  ncols = other.ncols;
  n = other.n;
  order = other.order;
  data = new T[n];
  std::copy(other.data, other.data + n, data);
  return *this;
}

template <class T>
Matrix<T>::Matrix(Matrix<T> &&other) noexcept
    : nrows(other.nrows), ncols(other.ncols),
      n(other.n), order(other.order), data(other.data) {
  other.data = nullptr;
}

template <class T>
Matrix<T> &Matrix<T>::operator=(Matrix<T> &&other) noexcept {
  if (this == &other)
    return *this;
  delete[] data;
  nrows = other.nrows;
  ncols = other.ncols;
  n = other.n;
  order = other.order;
  data = other.data;
  other.data = nullptr;
  return *this;
}

template <class T>
Matrix<T>::~Matrix() {
  delete[] data;
}

template <class T>
T &Matrix<T>::operator()(std::size_t row, std::size_t col) {
  if (!in_bounds(row, col))
    throw std::out_of_range("Matrix indices out of bounds");
  return data[index(row, col)];
}

template <class T>
const T &Matrix<T>::operator()(std::size_t row, std::size_t col) const {
  if (!in_bounds(row, col))
    throw std::out_of_range("Matrix indices out of bounds");
  return data[index(row, col)];
}

template <class T>
size_t Matrix<T>::size() const {
  return n;
}

template <class T>
std::size_t Matrix<T>::rows() const { return nrows; }

template <class T>
std::size_t Matrix<T>::cols() const { return ncols; }

template <class T>
std::size_t Matrix<T>::index(std::size_t row, std::size_t col) const {
  switch (order) {
  case ColMajorOrder:
    return row + col * nrows;  // column-major: (row + col * nrows)
  case RowMajorOrder:
    return row * ncols + col;  // row-major: (row * ncols + col)
  }
  // unreachable fallback:
  return row * ncols + col;
}

template <class T>
bool Matrix<T>::in_bounds(std::size_t row, std::size_t col) const {
  return row < nrows && col < ncols;
}

// Explicit instantiations
template class Matrix<int>;
template class Matrix<double>;

} // namespace beam
