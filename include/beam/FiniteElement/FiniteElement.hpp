#include <beam/FiniteElement/Geometry.hpp>

namespace beam {
class BasisType {
public:
  enum {
    Invalid = -1,
    GaussLegendre = 0,
    GaussLobatto = 1,
    Positive = 2,
    OpenUniform = 3,
    ClosedUniform = 4,
    NumBasisTypes = 5
  };
};

class FunctionSpace {
public:
  enum {
    Pk, ///< Polynomials of order k
  };
};

class FiniteElement {
protected:
  int dim;  ///< Dimension of reference space
  int vdim; ///< Vector dimension of vector-valued basis functions
  int cdim; ///< Dimension of curl for vector-valued basis functions
  Geometry::Type geom_type; ///< Geometry::Type of the reference element
  int func_space, range_type, map_type, deriv_type, deriv_range_type,
      deriv_map_type;
  mutable int dof, ///< Number of degrees of freedom
      order;       ///< Order/dergree of the shape functions
  // mutable int orders[Geometry::MaxDim];
  // IntegrationRule Nodes;

  enum RangeType { UNKNOWN_RANGE_TYPE = -1, SCALAR, VECTOR };
  enum MapType { UNKNOWN_MAP_TYPE = -1, VALUE, INTEGRAL, H_DIV, H_CURL };

  FiniteElement(int D, Geometry::Type G, int Do, int O, int F = FunctionSpace::Pk);
};
} // namespace beam