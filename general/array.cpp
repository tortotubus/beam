#include "../config/config.hpp"
#include "array.hpp"

namespace beam {

template <class T> T Array<T>::Max() const {
  // MFEM_ASSERT(size > 0, "Array is empty with size " << size);

  T max = operator[](0);
  for (int i = 1; i < size; i++) {
    if (max < operator[](i)) {
      max = operator[](i);
    }
  }

  return max;
}

template <class T> T Array<T>::Min() const {
  // MFEM_ASSERT(size > 0, "Array is empty with size " << size);

  T min = operator[](0);
  for (int i = 1; i < size; i++) {
    if (operator[](i) < min) {
      min = operator[](i);
    }
  }

  return min;
}

// Partial Sum
template <class T> void Array<T>::PartialSum() {
  T sum = static_cast<T>(0);
  for (int i = 0; i < size; i++) {
    sum += operator[](i);
    operator[](i) = sum;
  }
}

// Sum
template <class T> T Array<T>::Sum() const {
  T sum = static_cast<T>(0);
  for (int i = 0; i < size; i++) {
    sum += operator[](i);
  }

  return sum;
}

template <class T> int Array<T>::IsSorted() const {
  T val_prev = operator[](0), val;
  for (int i = 1; i < size; i++) {
    val = operator[](i);
    if (val < val_prev) {
      return 0;
    }
    val_prev = val;
  }

  return 1;
}

template class Array<char>;
template class Array<int>;
template class Array<long long>;
template class Array<real_t>;

} // namespace beam