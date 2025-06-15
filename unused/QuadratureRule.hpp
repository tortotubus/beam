#pragma once

#include <vector>
#include <cmath>
#include <stdexcept>

namespace beam {
class QuadratureRule {
public:
    QuadratureRule(size_t order) {
        switch(order) {
            case 1:
                // Midpoint rule
                xi = {0.0};
                w = {2.0};
                break;
            case 2:
                // Gauss-Legendre 2-point
                xi = {-1.0/std::sqrt(3.0), 1.0/std::sqrt(3.0)};
                w = {1.0, 1.0};
                break;
            case 3:
                // Gauss-Legendre 3-point
                xi = {-std::sqrt(0.6), 0.0, std::sqrt(0.6)};
                w = {5.0/9.0, 8.0/9.0, 5.0/9.0};
                break;
            default:
                throw std::runtime_error("Unsupported quadrature order");
        }
    }

    const std::vector<double>& points() const { return xi; }
    const std::vector<double>& weights() const { return w; }

private:
    std::vector<double> xi;
    std::vector<double> w;
};
} // namespace beam