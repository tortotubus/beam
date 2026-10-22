#pragma once
#include <functional>
#include <optional>
#include "elff/general/random.hpp"
#include "elff/c/models/ibm/IBCoupling.h"
#include "elff/contact/BroadPhaseLinkedCell.hpp"
#include "elff/contact/ContactGeometryProvider.hpp"
#include "elff/models/runtime/ModelRuntime.hpp"

namespace ELFF::Models {
class IBCoupling;

/** @brief Outcome of a setup-only geometry query or registration. */
enum class InsertionStatus {
  Accepted,
  ClearanceViolation,
  InvalidInput,
  UnsupportedGeometry,
  OutOfDomain,
  InvalidLifecycle,
  CollectiveMismatch,
  PreparationFailed
};

/** @brief Query/adoption result; model_id is assigned only by successful
 * adoption. Conflict indices identify the first detected violation, not a
 * global minimum.
 */
struct InsertionResult {
  InsertionStatus status = InsertionStatus::Accepted;
  ModelId model_id = invalid_model_id;
  ModelId conflicting_model = invalid_model_id;
  size_t candidate_segment = 0;
  size_t existing_segment = 0;
  bool accepted() const { return status == InsertionStatus::Accepted; }
};

/** @brief Fixed nonperiodic dense-grid domain for insertion centerlines. */
struct InsertionGridConfig {
  std::array<real_t, 3> origin;
  real_t cell_width;
  std::array<size_t, 3> dimensions;
};

/** @brief Setup service checking polylines before runtime adoption.
 * Owns fixed geometry snapshots; finalize model placement before registration.
 * No stepping, contact forces, self-clearance or checkpoint payload is
 * provided.
 */
class ModelInsertion : public RuntimeComponent {
public:
  using RuntimeComponent::RuntimeComponent;
  using CandidateFactory = std::function<std::unique_ptr<Model>(const CounterRandom&, std::uint64_t)>;
  /** @brief Configure locally once; MPI callers supply identical grids on all ranks.
   * Optional exact upper bounds restrict geometry inside a rounded grid; an
   * inactive axis may have zero extent. Search bounds still use the full grid.
   */
  InsertionStatus configure(const InsertionGridConfig &config,
      std::optional<std::array<real_t, 3>> exact_upper = std::nullopt);
  /** @brief Local check against snapshots; equality violates clearance.
   * This does not perform MPI agreement or adopt a model. */
  InsertionResult
  check_clearance(const Contact::ContactGeometryProvider &provider,
                  real_t clearance) const;
  /** @brief Serial-only snapshot of a local model without testing prior obstacles.
   */
  InsertionStatus register_model(ModelId id);
  /** @brief Collectively snapshot the owner's geometry on every rank.
   * No clearance check is performed. All ranks must supply the same ID.
   */
  InsertionStatus register_model_collective(ModelId id);
  /** @brief Serial-only adoption; rejection destroys the candidate without an ID.
   */
  InsertionResult try_adopt(ModelType type, std::unique_ptr<Model> candidate,
                            real_t clearance);
  /** @brief Adopt and register IB coupling and insertion geometry together. */
  InsertionResult try_adopt_ib(IBCoupling &coupling, ModelType type,
                               std::unique_ptr<Model> candidate,
                               real_t clearance,
                               ib_model_info_t info = {0, 0, -1, 0});
  /** @brief Collectively construct, check and adopt identical candidates.
   * Every communicator rank must call in the same order, with identical input
   * and a factory producing the same model. Factory exceptions become collective
   * PreparationFailed results. Only owner keeps the model; all keep snapshots.
   * The factory receives the component's fixed seed-1 placement stream and the
   * attempt index. Each invoked factory consumes one attempt, even if rejected
   * or throwing. Invalid collective setup consumes none. No external RNG setup
   * is needed. Also works serially (owner zero), including MPI-disabled builds.
   */
  InsertionResult try_adopt_replicated(ModelType type,
      const CandidateFactory& factory, real_t clearance,
      int owner, IBCoupling* coupling = nullptr,
      ib_model_info_t info = {0, 0, -1, 0});
  /** @brief End insertion; runtime sealing also invokes this hook. */
  int seal() override;
  size_t size() const { return snapshots_.size(); }
  std::uint64_t attempts() const { return next_attempt_; }

private:
  struct Snapshot {
    ModelId id;
    std::vector<std::array<real_t, 3>> points;
  };
  InsertionStatus agree_status(InsertionStatus) const;
  bool available() const;
  InsertionStatus validate(Contact::PolylineView geometry) const;
  InsertionResult query(Contact::PolylineView geometry, real_t clearance) const;
  InsertionResult adopt(ModelType type, std::unique_ptr<Model> candidate,
                        real_t clearance, IBCoupling *coupling,
                        ib_model_info_t info);
  const CounterRandom random_{1, 0x3333333333333333ULL};
  std::uint64_t next_attempt_ = 0;
  InsertionGridConfig config_{};
  std::array<real_t, 3> upper_{};
  std::array<real_t, 3> geometry_upper_{};
  bool configured_ = false;
  bool sealed_ = false;
  std::vector<Snapshot> snapshots_;
};

} // namespace ELFF::Models
