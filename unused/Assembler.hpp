#pragma once

#include <beam/FiniteElement/DoFHandler.hpp>
#include <beam/FiniteElement/QuadratureRule.hpp>
#include <functional>
#include <Eigen/Dense>

namespace beam {

template<std::size_t Dim, std::size_t NNodes>
class Assembler {
public:
    using BilinearForm = std::function<double(
        const std::vector<double>&,        // Values of test function
        const std::vector<double>&,        // Values of trial function
        const std::vector<double>&,        // Derivatives of test function
        const std::vector<double>&,        // Derivatives of trial function
        const std::vector<double>&,        // Second derivatives of test function
        const std::vector<double>&,        // Second derivatives of trial function
        double,                            // Quadrature point
        const std::array<double, Dim>&     // Physical point
    )>;

    using LinearForm = std::function<double(
        const std::vector<double>&,        // Values of test function
        const std::vector<double>&,        // Derivatives of test function
        const std::vector<double>&,        // Second derivatives of test function
        double,                            // Quadrature point
        const std::array<double, Dim>&     // Physical point
    )>;

    Assembler(const DoFHandler<Dim, NNodes>& dof_handler, 
              const QuadratureRule& quad_rule)
        : dof_handler_(dof_handler), quad_rule_(quad_rule) {}

    void assemble_matrix(const std::string& test_field,
                        const std::string& trial_field,
                        const BilinearForm& form,
                        Eigen::MatrixXd& matrix) {
        const auto& mesh = dof_handler_.get_mesh();
        const auto test_space = dof_handler_.get_field_space(test_field);
        const auto trial_space = dof_handler_.get_field_space(trial_field);

        for (size_t e = 0; e < mesh.size_elements(); ++e) {
            // Get global DOF indices for this element
            auto test_dofs = dof_handler_.element_dof_indices(e, test_field);
            auto trial_dofs = dof_handler_.element_dof_indices(e, trial_field);

            // For each quadrature point
            for (size_t q = 0; q < quad_rule_.points().size(); ++q) {
                double xi = quad_rule_.points()[q];
                double w = quad_rule_.weights()[q];

                // Evaluate basis functions
                std::vector<double> test_vals, trial_vals;
                std::vector<double> test_derivs, trial_derivs;
                std::vector<double> test_second_derivs, trial_second_derivs;

                test_space->evaluateAt(xi, test_vals);
                trial_space->evaluateAt(xi, trial_vals);
                test_space->evaluateDerivativeAt(xi, test_derivs);
                trial_space->evaluateDerivativeAt(xi, trial_derivs);
                test_space->evaluateSecondDerivativeAt(xi, test_second_derivs);
                trial_space->evaluateSecondDerivativeAt(xi, trial_second_derivs);

                // Get physical coordinates and Jacobian
                auto [x, J] = compute_physical_mapping(e, xi);

                // Scale derivatives by Jacobian
                for (auto& d : test_derivs) d /= J;
                for (auto& d : trial_derivs) d /= J;
                for (auto& d : test_second_derivs) d /= (J * J);
                for (auto& d : trial_second_derivs) d /= (J * J);

                // Evaluate bilinear form
                double val = form(test_vals, trial_vals, 
                                test_derivs, trial_derivs,
                                test_second_derivs, trial_second_derivs,
                                xi, x);

                // Add to global matrix
                for (size_t i = 0; i < test_dofs.size(); ++i) {
                    for (size_t j = 0; j < trial_dofs.size(); ++j) {
                        matrix(test_dofs[i], trial_dofs[j]) += val * w * J;
                    }
                }
            }
        }
    }

    void assemble_vector(const std::string& test_field,
                        const LinearForm& form,
                        Eigen::VectorXd& vector) {
        const auto& mesh = dof_handler_.get_mesh();
        const auto test_space = dof_handler_.get_field_space(test_field);

        for (size_t e = 0; e < mesh.size_elements(); ++e) {
            // Get global DOF indices for this element
            auto test_dofs = dof_handler_.element_dof_indices(e, test_field);

            // For each quadrature point
            for (size_t q = 0; q < quad_rule_.points().size(); ++q) {
                double xi = quad_rule_.points()[q];
                double w = quad_rule_.weights()[q];

                // Evaluate basis functions
                std::vector<double> test_vals;
                std::vector<double> test_derivs;
                std::vector<double> test_second_derivs;

                test_space->evaluateAt(xi, test_vals);
                test_space->evaluateDerivativeAt(xi, test_derivs);
                test_space->evaluateSecondDerivativeAt(xi, test_second_derivs);

                // Get physical coordinates and Jacobian
                auto [x, J] = compute_physical_mapping(e, xi);

                // Scale derivatives by Jacobian
                for (auto& d : test_derivs) d /= J;
                for (auto& d : test_second_derivs) d /= (J * J);

                // Evaluate linear form
                double val = form(test_vals, test_derivs, test_second_derivs, xi, x);

                // Add to global vector
                for (size_t i = 0; i < test_dofs.size(); ++i) {
                    vector(test_dofs[i]) += val * w * J;
                }
            }
        }
    }

private:
    // Returns physical point and Jacobian of the mapping
    std::pair<std::array<double, Dim>, double> 
    compute_physical_mapping(size_t elem_idx, double xi) const {
        const auto& mesh = dof_handler_.get_mesh();
        const auto& element = mesh.get_element(elem_idx);
        std::array<double, Dim> x{};
        
        // Linear interpolation between nodes
        double N1 = 0.5 * (1 - xi);
        double N2 = 0.5 * (1 + xi);
        double dN1 = -0.5;
        double dN2 = 0.5;

        const auto& x1 = mesh.get_node(element[0]);
        const auto& x2 = mesh.get_node(element[1]);

        // Compute physical point and Jacobian
        std::array<double, Dim> dx{};
        for (size_t d = 0; d < Dim; ++d) {
            x[d] = N1 * x1[d] + N2 * x2[d];
            dx[d] = dN1 * x1[d] + dN2 * x2[d];
        }

        // Compute Jacobian (length of tangent vector for curve in Dim dimensions)
        double J = 0.0;
        for (size_t d = 0; d < Dim; ++d) {
            J += dx[d] * dx[d];
        }
        J = std::sqrt(J);

        return {x, J};
    }

    const DoFHandler<Dim, NNodes>& dof_handler_;
    const QuadratureRule& quad_rule_;
};

} // namespace beam
