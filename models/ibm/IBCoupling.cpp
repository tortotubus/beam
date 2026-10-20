#include "elff/models/ibm/IBCoupling.hpp"

namespace ELFF::Models {
IBCoupling::Registration IBCoupling::describe(Model* model, ModelId id,
                                               ib_model_info_t info) {
  auto* ib = dynamic_cast<IBModel*>(model);
  auto* force = dynamic_cast<IBForceCoupled*>(ib);
  auto* velocity = dynamic_cast<IBVelocityCoupled*>(ib);
  if ((force == nullptr) == (velocity == nullptr))
    return {id, IB_MODEL_VIEW_INVALID, nullptr, info};
  return {id, force ? IB_MODEL_VIEW_FORCE_COUPLED : IB_MODEL_VIEW_VELOCITY_COUPLED,
          force ? static_cast<void*>(force) : static_cast<void*>(velocity), info};
}

int IBCoupling::register_model(ModelId id, ib_model_info_t info,
                               ib_model_view_kind_t remote_kind) {
  if (!runtime_.mutable_registration() || !valid_info(info)) return -1;
  const auto* entry = runtime_.entry(id);
  if (!entry) return -1;
  for (const auto& registered : registrations_)
    if (registered.id == id) return -1;
  Registration registration = describe(entry->model.get(), id, info);
  if (!entry->model) {
    if (remote_kind != IB_MODEL_VIEW_FORCE_COUPLED &&
        remote_kind != IB_MODEL_VIEW_VELOCITY_COUPLED) return -1;
    registration.kind = remote_kind;
  }
  if (registration.kind == IB_MODEL_VIEW_INVALID) return -1;
  registration.info.owner_rank = entry->owner_pid;
  registrations_.push_back(registration);
  return 0;
}

ModelId IBCoupling::add_remote_model(ModelType type, ib_model_view_kind_t kind,
                                    ib_model_info_t info) {
  if (!runtime_.mutable_registration() || !valid_info(info) ||
      (kind != IB_MODEL_VIEW_FORCE_COUPLED && kind != IB_MODEL_VIEW_VELOCITY_COUPLED))
    return invalid_model_id;
  registrations_.reserve(registrations_.size() + 1);
  if (!runtime_.reserve_model()) return invalid_model_id;
  const auto id = runtime_.add_remote_model(type, info.owner_rank);
  if (id == invalid_model_id) return id;
  registrations_.push_back({id, kind, nullptr, info});
  return id;
}

int IBCoupling::seal() {
  if (runtime_.sealed()) return 0;
  // Reserve before sorting so allocation failure leaves existing views alone.
  views_.reserve(registrations_.size());
  std::sort(registrations_.begin(), registrations_.end(),
            [](const auto& a, const auto& b) { return a.id < b.id; });
  views_.clear();
  for (auto& registration : registrations_) {
    registration.info.owner_rank = runtime_.entry(registration.id)->owner_pid;
    views_.push_back({registration.kind, registration.handle, &registration.info});
  }
  return 0;
}

int IBCoupling::get_model_views(ib_model_view_array_t& views) {
  views = {};
  if (runtime_.seal() != 0) return -1;
  views = {views_.size(), views_.empty() ? nullptr : views_.data()};
  return 0;
}
} // namespace ELFF::Models
