#include "fe_base.hpp"
#include "coefficient.hpp"
#include "../mesh/element.hpp"

namespace beam {
  namespace fem {
class BilinearFormIntegrator {
public:
  virtual ~BilinearFormIntegrator() = default;

  // virtual void AssembleElementMatrix(const FiniteElement &fe, const Element &el,
  //                                    const Coefficient &coeff,
  //                                    DenseMatrix &Ae) const = 0;
};
  }
} // namespace beam