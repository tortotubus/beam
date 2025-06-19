#pragma once

#include "../config/config.hpp"

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <initializer_list>
#include <iostream>
#include <type_traits>

namespace beam {

template <class T> class Array;

template <class T> void Swap(Array<T> &, Array<T> &);

/**
   Abstract data type Array.

   Array<T> is an automatically increasing array containing elements of the
   generic type T, which must be a trivial type, see `std::is_trivial`. The
   allocated size may be larger then the logical size of the array. The elements
   can be accessed by the [] operator, the range is 0 to size-1.
*/

template <class T> class Array {
protected:
  /// Pointer to data
  T *data;
  /// Size of the array
  int size, capacity;

  inline void GrowSize(int minsize);

  static_assert(std::is_trivial<T>::value, "type T must be trivial");

public:
  friend void Swap<T>(Array<T> &, Array<T> &);

  /// Creates an empty array
  inline Array() : size(0) {}

  /// Creates array of @a asize elements
  explicit inline Array(int asize) : size(asize) {
    data = new T[asize](); // Allocate and zero-initialize
  }

  /** @brief Creates array using an externally allocated host pointer @a data_
      to @a asize elements. If @a own_data is true, the array takes ownership
      of the pointer.

      When @a own_data is true, the pointer @a data_ must be allocated with
      MemoryType given by MemoryManager::GetHostMemoryType(). */
  inline Array(T *data_, int asize, bool own_data = false) {
    if (own_data) {
      data = data_;
    } else {
      data = new T[asize];
      std::copy(data_, data_ + asize, data);
    }
    size = asize;
  }

  /// Copy constructor: deep copy from @a src
  /** This method supports source arrays using any MemoryType. */
  inline Array(const Array &src) : size(src.size) {
    data = new T[size];
    std::copy(src.data, src.data + size, data);
  }

  /// Copy constructor (deep copy) from 'src', an Array of convertible type.
  template <typename CT> inline Array(const Array<CT> &src);

  /// Construct an Array from a C-style array of static length
  template <typename CT, int N> explicit inline Array(const CT (&values)[N]);

  /// Construct an Array from a braced initializer list of convertible type
  template <typename CT, typename std::enable_if<std::is_convertible<CT, T>::value, bool>::type = true>
  explicit inline Array(std::initializer_list<CT> values);

  /// Move constructor ("steals" data from 'src')
  inline Array(Array<T> &&src) : Array() { Swap(src, *this); }

  /// Destructor
  inline ~Array() { delete[] data; }

  /// Assignment operator: deep copy from 'src'.
  Array<T> &operator=(const Array<T> &src) {
    src.Copy(*this);
    return *this;
  }

  /// Assignment operator (deep copy) from @a src, an Array of convertible type.
  template <typename CT> inline Array &operator=(const Array<CT> &src);

  /// Return the data as 'T *'
  inline operator T *() { return data; }

  /// Return the data as 'const T *'
  inline operator const T *() const { return data; }

  /// Returns the data
  inline T *GetData() { return data; }

  /// Returns the data
  inline const T *GetData() const { return data; }

  /// Changes the ownership of the data
  inline void StealData(T **p) {
    // TODO
  }

  /// NULL-ifies the data
  inline void LoseData() {
    // TODO
  }

  /// Make the Array own the data
  void MakeDataOwner() const {
    // TODO
  }

  /// Return the logical size of the array.
  inline int Size() const { return size; }

  /// Change the logical size of the array, keep existing entries.
  inline void SetSize(int nsize);

  /// Same as SetSize(int) plus initialize new entries with 'initval'.
  inline void SetSize(int nsize, const T &initval);

  /** Maximum number of entries the array can store without allocating more
      memory. */
  inline int Capacity() const { return capacity; }

  /// Ensures that the allocated size is at least the given size.
  inline void Reserve(int capacity) {
    if (capacity > Capacity()) {
      GrowSize(capacity);
    }
  }

  /// Reference access to the ith element.
  inline T &operator[](int i);

  /// Const reference access to the ith element.
  inline const T &operator[](int i) const;

  /// Append element 'el' to array, resize if necessary.
  inline int Append(const T &el);

  /// Append another array to this array, resize if necessary.
  inline int Append(const T *els, int nels);

  /// Append another array to this array, resize if necessary.
  inline int Append(const Array<T> &els) { return Append(els, els.Size()); }

  /// Prepend an 'el' to the array, resize if necessary.
  inline int Prepend(const T &el);

  /// Return the last element in the array.
  inline T &Last();

  /// Return the last element in the array.
  inline const T &Last() const;

  /// Append element when it is not yet in the array, return index.
  inline int Union(const T &el);

  /// Return the first index where 'el' is found; return -1 if not found.
  inline int Find(const T &el) const;

  /// Do bisection search for 'el' in a sorted array; return -1 if not found.
  inline int FindSorted(const T &el) const;

  /// Delete the last entry of the array.
  inline void DeleteLast() {
    if (size > 0) {
      size--;
    }
  }

  /// Delete the first entry with value == 'el'.
  inline void DeleteFirst(const T &el);

  /// Delete the whole array.
  inline void DeleteAll();

  /// Reduces the capacity of the array to exactly match the current size.
  inline void ShrinkToFit();

  /// Create a copy of the internal array to the provided @a copy.
  inline void Copy(Array &copy) const;

  /// Make this Array a reference to a pointer.
  /** When @a own_data is true, the pointer @a data_ must be allocated with
      MemoryType given by MemoryManager::GetHostMemoryType(). */
  inline void MakeRef(T *data_, int size_, bool own_data = false);

  /// Make this Array a reference to 'master'.
  inline void MakeRef(const Array &master);

  /**
   * @brief Permute the array using the provided indices. Sorts the indices
   * variable in the process, thereby destroying the permutation. The rvalue
   * reference is to be used when this destruction is allowed, whilst the const
   * reference preserves at the cost of duplication.
   *
   * @param indices The indices of the ordering. data[i] = data[indices[i]].
   */
  template <typename I> inline void Permute(I &&indices);
  template <typename I> inline void Permute(const I &indices) {
    Permute(I(indices));
  }

  /// Copy sub array starting from @a offset out to the provided @a sa.
  inline void GetSubArray(int offset, int sa_size, Array<T> &sa) const;

  /** @brief Find the maximal element in the array, using the comparison
      operator `<` for class T. */
  T Max() const;

  /** @brief Find the minimal element in the array, using the comparison
      operator `<` for class T. */
  T Min() const;

  /// Sorts the array in ascending order. This requires operator< to be defined
  /// for T.
  void Sort() { std::sort((T *)data, data + size); }

  /// Sorts the array in ascending order using the supplied comparison function
  /// object.
  template <class Compare> void Sort(Compare cmp) {
    std::sort((T *)data, data + size, cmp);
  }

  /** @brief Removes duplicities from a sorted array. This requires
      operator== to be defined for T. */
  void Unique() {
    T *end = std::unique((T *)data, data + size);
    SetSize((int)(end - data));
  }

  /// Return 1 if the array is sorted from lowest to highest.  Otherwise return
  /// 0.
  int IsSorted() const;

  /// Does the Array have Size zero.
  bool IsEmpty() const { return Size() == 0; }

  /// Fill the entries of the array with the cumulative sum of the entries.
  void PartialSum();

  /// Return the sum of all the array entries using the '+'' operator for class
  /// 'T'.
  T Sum() const;

  /// Set all entries of the array to the provided constant.
  inline void operator=(const T &a);

  /// Copy data from a pointer. 'Size()' elements are copied.
  inline void Assign(const T *);

  /// STL-like copyTo @a dest from begin to end.
  template <typename U> inline void CopyTo(U *dest) {
    std::copy(begin(), end(), dest);
  }

  /** @brief Copy from @a src into this array.  Copies enough entries to
      fill the Capacity size of this array.  Careful this does not update
      the Size to match this Capacity after this.*/
  template <typename U> inline void CopyFrom(const U *src) {
    std::memcpy(begin(), src, MemoryUsage());
  }

  /// STL-like begin.  Returns pointer to the first element of the array.
  inline T *begin() { return data; }

  /// STL-like end.  Returns pointer after the last element of the array.
  inline T *end() { return data + size; }

  /// STL-like begin.  Returns const pointer to the first element of the array.
  inline const T *begin() const { return data; }

  /// STL-like end.  Returns const pointer after the last element of the array.
  inline const T *end() const { return data + size; }

  /// Returns the number of bytes allocated for the array including any reserve.
  std::size_t MemoryUsage() const { return Capacity() * sizeof(T); }
};

template <class T>
inline bool operator==(const Array<T> &LHS, const Array<T> &RHS) {
  if (LHS.Size() != RHS.Size()) {
    return false;
  }
  for (int i = 0; i < LHS.Size(); i++) {
    if (LHS[i] != RHS[i]) {
      return false;
    }
  }
  return true;
}

template <class T>
inline bool operator!=(const Array<T> &LHS, const Array<T> &RHS) {
  return !(LHS == RHS);
}

// /// Utility function similar to std::as_const in c++17.
// template <typename T> const T &AsConst(const T &a) { return a; }

/// inlines ///
template <class T> inline void Swap(T &a, T &b) {
  T c = a;
  a = b;
  b = c;
}

template <class T> inline void Swap(Array<T> &a, Array<T> &b) {
  Swap(a.data, b.data);
  Swap(a.size, b.size);
}

// template <class T> inline Array<T>::Array(const Array &src) :
// size(src.Size()) {
//   data = new T[size];
//   std::copy(src.data, src.data + size, data);
// }

template <typename T>
template <typename CT>
inline Array<T>::Array(const Array<CT> &src) : size(src.Size()) {
  data = new T[size];
  std::copy(src.GetData(), src.GetData() + size, data);
}

template <typename T>
template <typename CT, typename std::enable_if<
                           std::is_convertible<CT, T>::value, bool>::type>
inline Array<T>::Array(std::initializer_list<CT> values)
    : Array(values.size()) {
  std::copy(values.begin(), values.end(), begin());
}

template <typename T>
template <typename CT, int N>
inline Array<T>::Array(const CT (&values)[N]) : Array(N) {
  std::copy(values, values + N, begin());
}

template <class T> inline void Array<T>::GrowSize(int minsize) {
  if (minsize <= Capacity())
    return;
  int new_capacity = std::max(minsize, 2 * Capacity());
  T *new_data = new T[new_capacity];
  std::copy(data, data + size, new_data);
  delete[] data;
  data = new_data;
}

template <typename T> inline void Array<T>::ShrinkToFit() {
  if (Capacity() == size) { return; }
  T *temp = new T[size];
  std::copy(data, data + (size_t)size, temp);
  delete[] data;
  data = temp;
}

template <typename T>
template <typename I>
inline void Array<T>::Permute(I &&indices) {
  for (int i = 0; i < size; i++) {
    auto current = i;
    while (i != indices[current]) {
      auto next = indices[current];
      std::swap(data[current], data[next]);
      indices[current] = current;
      current = next;
    }
    indices[current] = current;
  }
}

template <typename T>
template <typename CT>
inline Array<T> &Array<T>::operator=(const Array<CT> &src) {
  SetSize(src.Size());
  for (int i = 0; i < size; i++) {
    (*this)[i] = T(src[i]);
  }
  return *this;
}

template <class T> inline void Array<T>::SetSize(int nsize) {
  // MFEM_ASSERT(nsize >= 0, "Size must be non-negative.  It is " << nsize);
  if (nsize > Capacity()) {
    GrowSize(nsize);
  }
  size = nsize;
}

template <class T> inline void Array<T>::SetSize(int nsize, const T &initval) {
  // MFEM_ASSERT(nsize >= 0, "Size must be non-negative.  It is " << nsize);
  if (nsize > size) {
    if (nsize > Capacity()) {
      GrowSize(nsize);
    }
    for (int i = size; i < nsize; i++) {
      data[i] = initval;
    }
  }
  size = nsize;
}

template <class T> inline T &Array<T>::operator[](int i) {
  // MFEM_ASSERT(i >= 0 && i < size, "Access element " << i << " of array, size
  // = " << size);
  return data[i];
}

template <class T> inline const T &Array<T>::operator[](int i) const {
  // MFEM_ASSERT(i >= 0 && i < size,  "Access element " << i << " of array, size
  // = " << size);
  return data[i];
}

template <class T> inline int Array<T>::Append(const T &el) {
  SetSize(size + 1);
  data[size - 1] = el;
  return size;
}

template <class T> inline int Array<T>::Append(const T *els, int nels) {
  const int old_size = size;

  SetSize(size + nels);
  for (int i = 0; i < nels; i++) {
    data[old_size + i] = els[i];
  }
  return size;
}

template <class T> inline int Array<T>::Prepend(const T &el) {
  SetSize(size + 1);
  for (int i = size - 1; i > 0; i--) {
    data[i] = data[i - 1];
  }
  data[0] = el;
  return size;
}

template <class T> inline T &Array<T>::Last() {
  // MFEM_ASSERT(size > 0, "Array size is zero: " << size);
  return data[size - 1];
}

template <class T> inline const T &Array<T>::Last() const {
  // MFEM_ASSERT(size > 0, "Array size is zero: " << size);
  return data[size - 1];
}

template <class T> inline int Array<T>::Union(const T &el) {
  int i = 0;
  while ((i < size) && (data[i] != el)) {
    i++;
  }
  if (i == size) {
    Append(el);
  }
  return i;
}

template <class T> inline int Array<T>::Find(const T &el) const {
  for (int i = 0; i < size; i++) {
    if (data[i] == el) {
      return i;
    }
  }
  return -1;
}

template <class T> inline int Array<T>::FindSorted(const T &el) const {
  const T *begin = data, *end = begin + size;
  const T *first = std::lower_bound(begin, end, el);
  if (first == end || !(*first == el)) {
    return -1;
  }
  return (int)(first - begin);
}

template <class T> inline void Array<T>::DeleteFirst(const T &el) {
  for (int i = 0; i < size; i++) {
    if (data[i] == el) {
      for (i++; i < size; i++) {
        data[i - 1] = data[i];
      }
      size--;
      return;
    }
  }
}

template <class T> inline void Array<T>::DeleteAll() {
  delete[] data;
  data = nullptr;
  size = 0;
}

template <typename T> inline void Array<T>::Copy(Array &copy) const {
  copy.SetSize(size);
  std::copy(data, data + size, copy.data);
}

template <class T>
inline void Array<T>::MakeRef(T *data_, int size_, bool own_data) {
  if (own_data) {
    delete[] data;
  }
  data = data_;
  size = size_;
}

template <class T> inline void Array<T>::MakeRef(const Array &master) {
  data = master.data;
  size = master.size;
}

template <class T>
inline void Array<T>::GetSubArray(int offset, int sa_size, Array<T> &sa) const {
  sa.SetSize(sa_size);
  for (int i = 0; i < sa_size; i++) {
    sa[i] = (*this)[offset + i];
  }
}

template <class T> inline void Array<T>::operator=(const T &a) {
  for (int i = 0; i < size; i++) {
    data[i] = a;
  }
}

template <class T> inline void Array<T>::Assign(const T *p) {
  std::copy(p, p + size, data);
}

} // namespace beam
