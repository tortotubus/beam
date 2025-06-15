#pragma once

#include <memory>
#include <vector>
#include "ShapeFunctionSet.hpp"

namespace beam {

class FunctionSpace {
public:
    virtual ~FunctionSpace() = default;
    virtual size_t numDofs() const = 0;
    virtual void evaluateAt(double xi, std::vector<double>& values) const = 0;
    virtual void evaluateDerivativeAt(double xi, std::vector<double>& derivatives) const = 0;
    virtual void evaluateSecondDerivativeAt(double xi, std::vector<double>& derivatives) const = 0;
};

class SingleSpace : public FunctionSpace {
public:
    SingleSpace(ShapeFunctionSet basis, size_t component_dim) 
        : basis_(basis), component_dim_(component_dim) {}
    
    size_t numDofs() const override {
        return basis_.num_basis_functions() * component_dim_;
    }

    void evaluateAt(double xi, std::vector<double>& values) const override {
        values.resize(numDofs());
        for (size_t i = 0; i < basis_.num_basis_functions(); ++i) {
            double val = basis_.evaluate(i, xi);
            for (size_t d = 0; d < component_dim_; ++d) {
                values[i * component_dim_ + d] = val;
            }
        }
    }

    void evaluateDerivativeAt(double xi, std::vector<double>& derivatives) const override {
        derivatives.resize(numDofs());
        for (size_t i = 0; i < basis_.num_basis_functions(); ++i) {
            double val = basis_.evaluate_derivative(i, xi);
            for (size_t d = 0; d < component_dim_; ++d) {
                derivatives[i * component_dim_ + d] = val;
            }
        }
    }

    void evaluateSecondDerivativeAt(double xi, std::vector<double>& derivatives) const override {
        derivatives.resize(numDofs());
        for (size_t i = 0; i < basis_.num_basis_functions(); ++i) {
            double val = basis_.evaluate_second_derivative(i, xi);
            for (size_t d = 0; d < component_dim_; ++d) {
                derivatives[i * component_dim_ + d] = val;
            }
        }
    }

private:
    ShapeFunctionSet basis_;
    size_t component_dim_;
};

class MixedSpace : public FunctionSpace {
public:
    void addSubspace(std::shared_ptr<FunctionSpace> space, size_t offset) {
        subspaces_.push_back({space, offset});
        total_dofs_ = std::max(total_dofs_, offset + space->numDofs());
    }
    
    size_t numDofs() const override {
        return total_dofs_;
    }

    void evaluateAt(double xi, std::vector<double>& values) const override {
        values.resize(total_dofs_);
        std::fill(values.begin(), values.end(), 0.0);
        std::vector<double> subvalues;
        
        for (const auto& [space, offset] : subspaces_) {
            space->evaluateAt(xi, subvalues);
            for (size_t i = 0; i < subvalues.size(); ++i) {
                values[offset + i] = subvalues[i];
            }
        }
    }

    void evaluateDerivativeAt(double xi, std::vector<double>& derivatives) const override {
        derivatives.resize(total_dofs_);
        std::fill(derivatives.begin(), derivatives.end(), 0.0);
        std::vector<double> subderivs;
        
        for (const auto& [space, offset] : subspaces_) {
            space->evaluateDerivativeAt(xi, subderivs);
            for (size_t i = 0; i < subderivs.size(); ++i) {
                derivatives[offset + i] = subderivs[i];
            }
        }
    }

    void evaluateSecondDerivativeAt(double xi, std::vector<double>& derivatives) const override {
        derivatives.resize(total_dofs_);
        std::fill(derivatives.begin(), derivatives.end(), 0.0);
        std::vector<double> subderivs;
        
        for (const auto& [space, offset] : subspaces_) {
            space->evaluateSecondDerivativeAt(xi, subderivs);
            for (size_t i = 0; i < subderivs.size(); ++i) {
                derivatives[offset + i] = subderivs[i];
            }
        }
    }

private:
    std::vector<std::pair<std::shared_ptr<FunctionSpace>, size_t>> subspaces_;
    size_t total_dofs_ = 0;
};

} // namespace beam
