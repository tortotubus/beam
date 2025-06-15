#pragma once

#include <beam/FiniteElement/Mesh.hpp>
#include <beam/FiniteElement/FieldVariable.hpp>
#include <beam/FiniteElement/FunctionSpace.hpp>
#include <map>
#include <string>
#include <memory>

namespace beam {

template <std::size_t DimRange, std::size_t NodesPerElement>
class DoFHandler {
public:
    DoFHandler(const Mesh<DimRange, NodesPerElement>& m) 
        : mesh_(m), total_dofs_(0) {}

    void register_field(const std::string& name, std::shared_ptr<FunctionSpace> space) {
        field_spaces_[name] = space;
        field_offsets_[name] = total_dofs_;
        total_dofs_ += space->numDofs();
    }

    size_t total_dofs() const { return total_dofs_; }

    std::vector<size_t> element_dof_indices(size_t elem_idx, const std::string& field_name) const {
        auto it_space = field_spaces_.find(field_name);
        auto it_offset = field_offsets_.find(field_name);
        
        if (it_space == field_spaces_.end() || it_offset == field_offsets_.end()) {
            throw std::runtime_error("Field not found: " + field_name);
        }

        const auto& element = mesh_.get_element(elem_idx);
        size_t dofs_per_node = it_space->second->numDofs() / NodesPerElement;
        std::vector<size_t> indices;
        indices.reserve(dofs_per_node * NodesPerElement);

        for (size_t node_idx = 0; node_idx < NodesPerElement; ++node_idx) {
            size_t node = element[node_idx];
            for (size_t dof = 0; dof < dofs_per_node; ++dof) {
                indices.push_back(it_offset->second + node * dofs_per_node + dof);
            }
        }

        return indices;
    }

    const Mesh<DimRange, NodesPerElement>& get_mesh() const { return mesh_; }

    std::shared_ptr<FunctionSpace> get_field_space(const std::string& field_name) const {
        auto it = field_spaces_.find(field_name);
        if (it == field_spaces_.end()) {
            throw std::runtime_error("Field not found: " + field_name);
        }
        return it->second;
    }

    size_t get_field_offset(const std::string& field_name) const {
        auto it = field_offsets_.find(field_name);
        if (it == field_offsets_.end()) {
            throw std::runtime_error("Field not found: " + field_name);
        }
        return it->second;
    }

private:
    const Mesh<DimRange, NodesPerElement>& mesh_;
    std::map<std::string, std::shared_ptr<FunctionSpace>> field_spaces_;
    std::map<std::string, size_t> field_offsets_;
    size_t total_dofs_;
};
} // namespace beam