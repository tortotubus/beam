#pragma once

#include <vector>
#include <array>
#include <stdexcept>
#include <cmath>

namespace beam {

class ShapeFunctionSet {
public:
    enum class Type { HermiteCubic, LagrangeLinear };

    ShapeFunctionSet(Type type) : basis_type_(type) {}

    size_t num_basis_functions() const {
        switch (basis_type_) {
            case Type::HermiteCubic:
                return 4;  // value and derivative at each endpoint
            case Type::LagrangeLinear:
                return 2;  // value at each endpoint
            default:
                throw std::runtime_error("Unknown basis type");
        }
    }

    double evaluate(size_t i, double xi) const {
        switch (basis_type_) {
            case Type::HermiteCubic:
                return evaluate_hermite(i, xi);
            case Type::LagrangeLinear:
                return evaluate_lagrange(i, xi);
            default:
                throw std::runtime_error("Unknown basis type");
        }
    }

    double evaluate_derivative(size_t i, double xi) const {
        switch (basis_type_) {
            case Type::HermiteCubic:
                return evaluate_hermite_derivative(i, xi);
            case Type::LagrangeLinear:
                return evaluate_lagrange_derivative(i, xi);
            default:
                throw std::runtime_error("Unknown basis type");
        }
    }

    double evaluate_second_derivative(size_t i, double xi) const {
        switch (basis_type_) {
            case Type::HermiteCubic:
                return evaluate_hermite_second_derivative(i, xi);
            case Type::LagrangeLinear:
                return 0.0;  // Linear functions have zero second derivative
            default:
                throw std::runtime_error("Unknown basis type");
        }
    }

    Type type() const { return basis_type_; }

private:
    Type basis_type_;

    double evaluate_hermite(size_t i, double xi) const {
        // Hermite basis functions on [-1, 1]
        switch(i) {
            case 0: return 0.25 * (2 - 3*xi + xi*xi*xi);  // H_00
            case 1: return 0.25 * (1 - xi - xi*xi + xi*xi*xi);  // H_10
            case 2: return 0.25 * (2 + 3*xi - xi*xi*xi);  // H_01
            case 3: return 0.25 * (-1 - xi + xi*xi + xi*xi*xi);  // H_11
            default: throw std::runtime_error("Invalid basis function index");
        }
    }

    double evaluate_lagrange(size_t i, double xi) const {
        // Linear Lagrange basis functions on [-1, 1]
        switch(i) {
            case 0: return 0.5 * (1 - xi);
            case 1: return 0.5 * (1 + xi);
            default: throw std::runtime_error("Invalid basis function index");
        }
    }

    double evaluate_hermite_derivative(size_t i, double xi) const {
        // Derivatives of Hermite basis functions
        switch(i) {
            case 0: return 0.25 * (-3 + 3*xi*xi);
            case 1: return 0.25 * (-1 - 2*xi + 3*xi*xi);
            case 2: return 0.25 * (3 - 3*xi*xi);
            case 3: return 0.25 * (-1 + 2*xi + 3*xi*xi);
            default: throw std::runtime_error("Invalid basis function index");
        }
    }

    double evaluate_lagrange_derivative(size_t i, double xi) const {
        // Derivatives of Linear Lagrange basis functions
        switch(i) {
            case 0: return -0.5;
            case 1: return 0.5;
            default: throw std::runtime_error("Invalid basis function index");
        }
    }

    double evaluate_hermite_second_derivative(size_t i, double xi) const {
        // Second derivatives of Hermite basis functions
        switch(i) {
            case 0: return 0.25 * (6*xi);
            case 1: return 0.25 * (-2 + 6*xi);
            case 2: return 0.25 * (-6*xi);
            case 3: return 0.25 * (2 + 6*xi);
            default: throw std::runtime_error("Invalid basis function index");
        }
    }
};

}