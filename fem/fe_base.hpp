#include "geom.hpp"

namespace beam {
  namespace fem {
class FiniteElement {
protected:
  int dim;  ///< Dimension of reference space
  int vdim; ///< Vector dimension of vector-valued basis functions
  int cdim; ///< Dimension of curl for vector-valued basis functions
  Geometry::Type geom_type; ///< Geometry::Type of the reference element
  int dof, order;

public:
  enum RangeType { UNKNOWN_RANGE_TYPE = -1, SCALAR, VECTOR };
  enum MapType { UNKNOWN_MAP_TYPE = -1, VALUE, INTEGRAL, H_DIV, H_CURL };
  enum DerivType { NONE, GRAD, DIV, CURL };

  //FiniteElement(int D, Geometry::Type G, int Do, int O, int F = FunctionSpace::Pk);
  FiniteElement(int O) : order(O) {};

  int GetDof() const { return dof; }
  int GetOrder() const { return order; }
};
  }
} // namespace beam