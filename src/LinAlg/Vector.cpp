#include <beam/LinAlg/Vector.hpp>
#include <algorithm>
#include <stdexcept>

namespace beam {

template <class T>
Vector<T>::Vector(size_t n) : n(n), data(new T[n]()) {
  // ...existing code...
}

template <class T>
Vector<T>::Vector(std::initializer_list<T> init_list)
    : n(init_list.size()), data(new T[init_list.size()]) {
  std::copy(init_list.begin(), init_list.end(), data);
}

template <class T>
Vector<T>::Vector(const Vector &other) : n(other.n) {
  data = new T[n];
  std::copy(other.data, other.data + n, data);
}

template <class T>
Vector<T> &Vector<T>::operator=(const Vector &other) {
  if (this == &other)
    return *this;
  delete[] data;
  n = other.n;
  data = new T[n];
  std::copy(other.data, other.data + n, data);
  return *this;
}

template <class T>
Vector<T>::Vector(Vector &&other) noexcept : n(other.n), data(other.data) {
  other.data = nullptr;
}

template <class T>
Vector<T> &Vector<T>::operator=(Vector &&other) noexcept {
  if (this == &other)
    return *this;
  delete[] data;
  n = other.n;
  data = other.data;
  other.data = nullptr;
  return *this;
}

template <class T>
Vector<T>::~Vector() {
  delete[] data;
}

template <class T>
T &Vector<T>::operator()(size_t i) {
  if (in_bounds(i) && data != nullptr) {
    return data[i];
  } else {
    throw std::out_of_range("Vector index out of bounds");
  }
}

template <class T>
const T &Vector<T>::operator()(size_t i) const {
  if (in_bounds(i) && data != nullptr) {
    return data[i];
  } else {
    throw std::out_of_range("Vector index out of bounds");
  }
}

template <class T>
size_t Vector<T>::size() const {
  return n;
}

template <class T>
inline bool Vector<T>::in_bounds(size_t i) const {
  return i < n;
}

// Explicit instantiations
template class Vector<int>;
template class Vector<double>;

}  // namespace beam

