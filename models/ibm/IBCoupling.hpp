#pragma once
#include "elff/models/ModelRuntime.hpp"
#include "elff/models/ibm/IBForceCoupled.hpp"
#include "elff/models/ibm/IBVelocityCoupled.hpp"
#include "elff/c/models/ibm/IBCoupling.h"

namespace ELFF::Models {
/**
 * @brief Runtime-owned IB service registering a subset of runtime models.
 *
 * Model IDs and ownership belong to ModelRuntime. This component owns only IB
 * metadata and borrowed coupling views. Views are ordered by ModelId and become
 * stable when the owning runtime is sealed.
 */
class IBCoupling : public RuntimeComponent {
public:
  explicit IBCoupling(ModelRuntime& runtime) : RuntimeComponent(runtime) {}

  /** @brief Register an existing local model, or a remote model with explicit kind. */
  int register_model(ModelId id, ib_model_info_t info = {0, 0, -1, 0},
                     ib_model_view_kind_t remote_kind = IB_MODEL_VIEW_INVALID);

  /** @brief Atomically adopt and register a local IB model; null indicates rejection. */
  template<class T>
  T* adopt_owned_model(ModelType type, std::unique_ptr<T> model,
                       ib_model_info_t info = {0, 0, -1, 0}) {
    static_assert(std::is_base_of_v<IBModel, T>);
    auto registration = describe(model.get(), runtime_.size(), info);
    if (!runtime_.mutable_registration() || type == ModelType::Invalid ||
        registration.kind == IB_MODEL_VIEW_INVALID || !valid_info(info)) return nullptr;
    registrations_.reserve(registrations_.size() + 1);
    if (!runtime_.reserve_model()) return nullptr;
    T* result = runtime_.adopt_owned_model(type, std::move(model));
    if (!result) return nullptr;
    registrations_.push_back(registration); // No allocations remain after adoption.
    return result;
  }

  /** @brief Atomically append a remote runtime entry and its IB metadata. */
  ModelId add_remote_model(ModelType type, ib_model_view_kind_t kind,
                           ib_model_info_t info);
  int get_model_views(ib_model_view_array_t& views);
  int seal() override;

private:
  struct Registration {
    ModelId id;
    ib_model_view_kind_t kind;
    void* handle;
    ib_model_info_t info;
  };
  static bool valid_info(const ib_model_info_t& info) {
    return (info.sparse == 0 || info.sparse == 1) && info.level >= -1 && info.node_count >= 0;
  }
  static Registration describe(Model* model, ModelId id, ib_model_info_t info);
  std::vector<Registration> registrations_;
  std::vector<ib_model_view_t> views_;
};
} // namespace ELFF::Models
