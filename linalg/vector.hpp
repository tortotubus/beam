#pragma once

#include "../general/array.hpp"
#include "../general/error.hpp"
#include "../general/globals.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <type_traits>
#include <initializer_list>


namespace beam
{

/// Vector data type.
class Vector
{
protected:
   real_t *data;
   int size, capacity;
   bool owns;
public:

   /** Default constructor for Vector. Sets size = 0, and calls Memory::Reset on
       data through Memory<double>'s default constructor. */
   Vector(): size(0), owns(true) {}

   /// Copy constructor. Allocates a new data array and copies the data.
   Vector(const Vector &);

   /// Move constructor. "Steals" data from its argument.
   Vector(Vector&& v);

   /// @brief Creates vector of size s.
   /// @warning Entries are not initialized to zero!
   explicit Vector(int s);

   /// Creates a vector referencing an array of doubles, owned by someone else.
   /** The pointer @a data_ can be NULL. The data array can be replaced later
       with SetData(). */
   //Vector(real_t *data_, int size_) { data.Wrap(data_, size_, false); size = size_; }

   /** @brief Create a Vector referencing a sub-vector of the Vector @a base
       starting at the given offset, @a base_offset, and size @a size_. */
   // Vector(Vector &base, int base_offset, int size_) : data(base.data, base_offset, size_), size(size_) { }

   /// Create a Vector of size @a size_ using MemoryType @a mt.
   // Vector(int size_, MemoryType mt) : data(size_, mt), size(size_) { }

   /** @brief Create a Vector of size @a size_ using host MemoryType @a h_mt and
       device MemoryType @a d_mt. */
   // Vector(int size_, MemoryType h_mt, MemoryType d_mt) : data(size_, h_mt, d_mt), size(size_) { }

   /// Create a vector from a statically sized C-style array of convertible type
   template <typename CT, int N>
   explicit Vector(const CT (&values)[N]) : Vector(N) { std::copy(values, values + N, begin()); }

   /// Create a vector using a braced initializer list
   template <typename CT, typename std::enable_if<std::is_convertible<CT,real_t>::value,bool>::type = true>
   explicit Vector(std::initializer_list<CT> values) : Vector(values.size()) { std::copy(values.begin(), values.end(), begin()); }

   /// Enable execution of Vector operations using the beam::Device.
   /** The default is to use Backend::CPU (serial execution on each MPI rank),
       regardless of the beam::Device configuration.

       When appropriate, BEAM functions and class methods will enable the use
       of the beam::Device for their Vector parameters.

       Some derived classes, e.g. GridFunction, enable the use of the
       beam::Device by default. */
   // virtual void UseDevice(bool use_dev) const { data.UseDevice(use_dev); }

   /// Return the device flag of the Memory object used by the Vector
   // virtual bool UseDevice() const { return data.UseDevice(); }

   /// Reads a vector from multiple files
   // void Load(std::istream ** in, int np, int * dim);

   /// Load a vector from an input stream.
   // void Load(std::istream &in, int Size);

   /// Load a vector from an input stream, reading the size from the stream.
   // void Load(std::istream &in) { int s; in >> s; Load(in, s); }

   /// @brief Resize the vector to size @a s.
   /** If the new size is less than or equal to Capacity() then the internal
       data array remains the same. Otherwise, the old array is deleted, if
       owned, and a new array of size @a s is allocated without copying the
       previous content of the Vector.
       @warning In the second case above (new size greater than current one),
       the vector will allocate new data array, even if it did not own the
       original data! Also, new entries are not initialized! */
   void SetSize(int s);

   /// Resize the vector to size @a s using MemoryType @a mt.
   // void SetSize(int s, MemoryType mt);

   /// Resize the vector to size @a s using the MemoryType of @a v.
   // void SetSize(int s, const Vector &v) { SetSize(s, v.GetMemory().GetMemoryType()); }

   /// Set the Vector data.
   /// @warning This method should be called only when OwnsData() is false.
   // void SetData(real_t *d) { data.Wrap(d, data.Capacity(), false); }

   /// Set the Vector data and size.
   /** The Vector does not assume ownership of the new data. The new size is
       also used as the new Capacity().
       @warning This method should be called only when OwnsData() is false.
       @sa NewDataAndSize(). */
   // void SetDataAndSize(real_t *d, int s) { data.Wrap(d, s, false); size = s; }

   /// Set the Vector data and size, deleting the old data, if owned.
   /** The Vector does not assume ownership of the new data. The new size is
       also used as the new Capacity().
       @sa SetDataAndSize(). */
   void NewDataAndSize(real_t *d, int s)
   { 
    // TODO
   }

   /// Reset the Vector to use the given external Memory @a mem and size @a s.
   /** If @a own_mem is false, the Vector will not own any of the pointers of
       @a mem.

       Note that when @a own_mem is true, the @a mem object can be destroyed
       immediately by the caller but `mem.Delete()` should NOT be called since
       the Vector object takes ownership of all pointers owned by @a mem.

       @sa NewDataAndSize(). */
   // inline void NewMemoryAndSize(const Memory<real_t> &mem, int s, bool own_mem);

   /// Reset the Vector to be a reference to a sub-vector of @a base.
   inline void MakeRef(Vector &base, int offset, int size);

   /** @brief Reset the Vector to be a reference to a sub-vector of @a base
       without changing its current size. */
   inline void MakeRef(Vector &base, int offset);

   /// Set the Vector data (host pointer) ownership flag.
   // void MakeDataOwner() const { data.SetHostPtrOwner(true); }

   /// Destroy a vector
   void Destroy();

   /// Returns the size of the vector.
   inline int Size() const { return size; }

   /// Return the size of the currently allocated data array.
   /** It is always true that Capacity() >= Size(). */
   inline int Capacity() const { return capacity; }

   /// Return a pointer to the beginning of the Vector data.
   /** @warning This method should be used with caution as it gives write access
       to the data of const-qualified Vector%s. */
   inline real_t *GetData() const { return const_cast<real_t*>((const real_t*)data); }

   /// STL-like begin.
   inline real_t *begin() { return data; }

   /// STL-like end.
   inline real_t *end() { return data + size; }

   /// STL-like begin (const version).
   inline const real_t *begin() const { return data; }

   /// STL-like end (const version).
   inline const real_t *end() const { return data + size; }

   /// Access Vector entries. Index i = 0 .. size-1.
   real_t &Elem(int i);

   /// Read only access to Vector entries. Index i = 0 .. size-1.
   const real_t &Elem(int i) const;

   /// Access Vector entries using () for 0-based indexing.
   /** @note If BEAM_DEBUG is enabled, bounds checking is performed. */
   inline real_t &operator()(int i);

   /// Read only access to Vector entries using () for 0-based indexing.
   /** @note If BEAM_DEBUG is enabled, bounds checking is performed. */
   inline const real_t &operator()(int i) const;

   /// Access Vector entries using [] for 0-based indexing.
   /** @note If BEAM_DEBUG is enabled, bounds checking is performed. */
   inline real_t &operator[](int i) { return (*this)(i); }

   /// Read only access to Vector entries using [] for 0-based indexing.
   /** @note If BEAM_DEBUG is enabled, bounds checking is performed. */
   inline const real_t &operator[](int i) const { return (*this)(i); }

   /// Dot product with a `double *` array.
   /// This function always executes on the CPU. A HostRead() will be called if
   /// required.
   /// To optionally execute on the device:
   /// Vector tmp(v, Size());
   /// res = (*this) * tmp;
   real_t operator*(const real_t *v) const;

   /// Return the inner-product.
   real_t operator*(const Vector &v) const;

   /// Copy Size() entries from @a v.
   Vector &operator=(const real_t *v);

   /// Copy assignment.
   /** @note Defining this method overwrites the implicitly defined copy
       assignment operator. */
   Vector &operator=(const Vector &v);

   /// Move assignment
   Vector &operator=(Vector&& v);

   /// Redefine '=' for vector = constant.
   Vector &operator=(real_t value);

   Vector &operator*=(real_t c);

   /// Component-wise scaling: (*this)(i) *= v(i)
   Vector &operator*=(const Vector &v);

   Vector &operator/=(real_t c);

   /// Component-wise division: (*this)(i) /= v(i)
   Vector &operator/=(const Vector &v);

   Vector &operator-=(real_t c);

   Vector &operator-=(const Vector &v);

   Vector &operator+=(real_t c);

   Vector &operator+=(const Vector &v);

   /// (*this) += a * Va
   Vector &Add(const real_t a, const Vector &Va);

   /// (*this) = a * x
   Vector &Set(const real_t a, const Vector &x);

   /// (*this)[i + offset] = v[i]
   void SetVector(const Vector &v, int offset);

   /// (*this)[i + offset] += v[i]
   void AddSubVector(const Vector &v, int offset);

   /// (*this) = -(*this)
   void Neg();

   /// (*this)(i) = 1.0 / (*this)(i)
   void Reciprocal();

   /// Swap the contents of two Vectors
   inline void Swap(Vector &other);

   /// Set v = v1 + v2.
   friend void add(const Vector &v1, const Vector &v2, Vector &v);

   /// Set v = v1 + alpha * v2.
   friend void add(const Vector &v1, real_t alpha, const Vector &v2, Vector &v);

   /// z = a * (x + y)
   friend void add(const real_t a, const Vector &x, const Vector &y, Vector &z);

   /// z = a * x + b * y
   friend void add(const real_t a, const Vector &x,
                   const real_t b, const Vector &y, Vector &z);

   /// Set v = v1 - v2.
   friend void subtract(const Vector &v1, const Vector &v2, Vector &v);

   /// z = a * (x - y)
   friend void subtract(const real_t a, const Vector &x,
                        const Vector &y, Vector &z);

   /// Computes cross product of this vector with another 3D vector.
   /// vout = this x vin.
   void cross3D(const Vector &vin, Vector &vout) const;

   /// v = median(v,lo,hi) entrywise.  Implementation assumes lo <= hi.
   void median(const Vector &lo, const Vector &hi);

   /// Extract entries listed in @a dofs to the output Vector @a elemvect.
   /** Negative dof values cause the -dof-1 position in @a elemvect to receive
       the -val in from this Vector. */
   void GetSubVector(const Array<int> &dofs, Vector &elemvect) const;

   /// Extract entries listed in @a dofs to the output array @a elem_data.
   /** Negative dof values cause the -dof-1 position in @a elem_data to receive
       the -val in from this Vector. */
   void GetSubVector(const Array<int> &dofs, real_t *elem_data) const;

   /// Set the entries listed in @a dofs to the given @a value.
   /** Negative dof values cause the -dof-1 position in this Vector to receive
       the -value. */
   void SetSubVector(const Array<int> &dofs, const real_t value);

   /// Set the entries listed in @a dofs to the given @a value (always on host).
   /** Negative dof values cause the -dof-1 position in this Vector to receive
       the -value.

       As opposed to SetSubVector(const Array<int>&, const real_t), this
       function will execute only on host, even if the vector or the @a dofs
       array have the device flag set. */
   void SetSubVectorHost(const Array<int> &dofs, const real_t value);

   /** @brief Set the entries listed in @a dofs to the values given in the @a
       elemvect Vector. Negative dof values cause the -dof-1 position in this
       Vector to receive the -val from @a elemvect. */
   void SetSubVector(const Array<int> &dofs, const Vector &elemvect);

   /** @brief Set the entries listed in @a dofs to the values given the @a ,
       elem_data array. Negative dof values cause the -dof-1 position in this
       Vector to receive the -val from @a elem_data. */
   void SetSubVector(const Array<int> &dofs, real_t *elem_data);

   /** @brief Add elements of the @a elemvect Vector to the entries listed in @a
       dofs. Negative dof values cause the -dof-1 position in this Vector to add
       the -val from @a elemvect. */
   void AddElementVector(const Array<int> & dofs, const Vector & elemvect);

   /** @brief Add elements of the @a elem_data array to the entries listed in @a
       dofs. Negative dof values cause the -dof-1 position in this Vector to add
       the -val from @a elem_data. */
   void AddElementVector(const Array<int> & dofs, real_t *elem_data);

   /** @brief Add @a times the elements of the @a elemvect Vector to the entries
       listed in @a dofs. Negative dof values cause the -dof-1 position in this
       Vector to add the -a*val from @a elemvect. */
   void AddElementVector(const Array<int> & dofs, const real_t a,
                         const Vector & elemvect);

   /// Set all vector entries NOT in the @a dofs Array to the given @a val.
   void SetSubVectorComplement(const Array<int> &dofs, const real_t val);

   /// Prints vector to stream out.
   // void Print(std::ostream &out = beam::out, int width = 8) const;

   /// Set random values in the vector.
   void Randomize(int seed = 0);
   /// Returns the l2 norm of the vector.
   real_t Norml2() const;
   /// Returns the l_infinity norm of the vector.
   real_t Normlinf() const;
   /// Returns the l_1 norm of the vector.
   real_t Norml1() const;
   /// Returns the l_p norm of the vector.
   real_t Normlp(real_t p) const;
   /// Returns the maximal element of the vector.
   real_t Max() const;
   /// Returns the minimal element of the vector.
   real_t Min() const;
   /// Return the sum of the vector entries
   real_t Sum() const;
   /// Compute the square of the Euclidean distance to another vector.
   inline real_t DistanceSquaredTo(const real_t *p) const;
   /// Compute the square of the Euclidean distance to another vector.
   inline real_t DistanceSquaredTo(const Vector &p) const;
   /// Compute the Euclidean distance to another vector.
   inline real_t DistanceTo(const real_t *p) const;
   /// Compute the Euclidean distance to another vector.
   inline real_t DistanceTo(const Vector &p) const;

   /// Destroys vector.
   virtual ~Vector();

};

// Inline methods

template <typename T>
inline T ZeroSubnormal(T val)
{
   return (std::fpclassify(val) == FP_SUBNORMAL) ? 0.0 : val;
}

inline bool IsFinite(const real_t &val)
{
   return std::isfinite(val);
}

inline int CheckFinite(const real_t *v, const int n)
{
   int bad = 0;
   for (int i = 0; i < n; i++)
   {
      if (!IsFinite(v[i])) { bad++; }
   }
   return bad;
}

inline Vector::Vector(int s)
{
  BEAM_ASSERT(s>=0,"Unexpected negative size.");
  size = s;
  data = new real_t[size]();
}

inline void Vector::SetSize(int s)
{
  if (s == size) {
    return;
  } 
  if (s <= capacity) {
    size = s;
    return;
  }
  size = s;
  real_t *new_data = new real_t[size];
  std::copy(data, data + size, new_data);
  delete[] data;
  data = new_data;
}

inline void Vector::MakeRef(Vector &base, int offset, int s)
{
  delete[] data;
  data = nullptr;

  size = s;
  data = base.data + offset;
}

inline void Vector::MakeRef(Vector &base, int offset)
{
   BEAM_ASSERT(offset >= 0 && offset <= base.size, "Offset [" << offset << "] is out of range [0," << base.size << "]");
  
  data = base.data + offset;
  size = base.size - offset;
}

inline void Vector::Destroy()
{
  size = 0;
  capacity = 0;
  if (owns) { delete[] data; }
  data = nullptr;
}

inline real_t &Vector::operator()(int i)
{
   BEAM_ASSERT(data && i >= 0 && i < size,
               "index [" << i << "] is out of range [0," << size << ")");

   return data[i];
}

inline const real_t &Vector::operator()(int i) const
{
   BEAM_ASSERT(data && i >= 0 && i < size,
               "index [" << i << "] is out of range [0," << size << ")");

   return data[i];
}

inline void Vector::Swap(Vector &other)
{
   beam::Swap(data, other.data);
   beam::Swap(size, other.size);
}

/// Specialization of the template function Swap<> for class Vector
template<> inline void Swap<Vector>(Vector &a, Vector &b)
{
   a.Swap(b);
}

inline Vector::~Vector()
{
  if (owns)
    delete[] data;
}

inline real_t DistanceSquared(const real_t *x, const real_t *y, const int n)
{
   real_t d = 0.0;

   for (int i = 0; i < n; i++)
   {
      d += (x[i]-y[i])*(x[i]-y[i]);
   }

   return d;
}

inline real_t Distance(const real_t *x, const real_t *y, const int n)
{
   return std::sqrt(DistanceSquared(x, y, n));
}

inline real_t Distance(const Vector &x, const Vector &y)
{
   return x.DistanceTo(y);
}

inline real_t Vector::DistanceSquaredTo(const real_t *p) const
{
   return DistanceSquared(data, p, size);
}

inline real_t Vector::DistanceSquaredTo(const Vector &p) const
{
   BEAM_ASSERT(p.Size() == Size(), "Incompatible vector sizes.");
   return DistanceSquared(data, p.data, size);
}

inline real_t Vector::DistanceTo(const real_t *p) const
{
   return Distance(data, p, size);
}

inline real_t Vector::DistanceTo(const Vector &p) const
{
   BEAM_ASSERT(p.Size() == Size(), "Incompatible vector sizes.");
   return Distance(data, p.data, size);
}

/// Returns the inner product of x and y
/** In parallel this computes the inner product of the local vectors,
    producing different results on each MPI rank.
*/
inline real_t InnerProduct(const Vector &x, const Vector &y)
{
   return x * y;
}


} // namespace beam
