#pragma once

#include <beam/FiniteElement/ShapeFunctionSet.hpp>

#include <string>
#include <vector>

namespace beam {
class FieldVariable {
public:
  FieldVariable(std::string name, size_t dofs_per_node, ShapeFunctionSet basis);

  size_t total_dofs() const;
  double &operator[](size_t global_index);
  const ShapeFunctionSet &shape_functions() const;
  const std::string name() const;

private:
  std::string name_;
  size_t dofs_per_node_;
  ShapeFunctionSet basis_;
  std::vector<double> coefficients_;
};
} // namespace beam