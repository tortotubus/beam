#pragma once

#include "../config/config.hpp"
#include "geom.hpp"
#include "fe_base.hpp"

#include <memory>

namespace beam {
  namespace fem {
class FiniteElementCollection {
  public:
    // Factory: Create FE for a given element type and order
    virtual std::unique_ptr<FiniteElement> Create(int element_type, int order) const = 0;
};
class CubicHermite1DCollection : public FiniteElementCollection {};

// class LinearFECollection : public FiniteElementCollection {
//   private: 
//    const PointFiniteElement PointFE;
//    const Linear1DFiniteElement SegmentFE;
//    const Linear2DFiniteElement TriangleFE;
//    const BiLinear2DFiniteElement QuadrilateralFE;
//    public:
//     LinearFECollection() : FiniteElementCollection(1) {}
// 
//     const FiniteElement * FiniteElementForGeometry(Geometry::Type GeomType) const override;
//     int DofForGeometry(Geometry::Type GeomType) const override;
// };
  }
} // namespace beam
