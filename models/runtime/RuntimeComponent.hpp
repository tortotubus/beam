#pragma once
namespace ELFF::Models {
class ModelRuntime;
/** @brief Optional runtime-owned service borrowing model identity and lifetime. */
class RuntimeComponent {
public:
  explicit RuntimeComponent(ModelRuntime& runtime) : runtime_(runtime) {}
  virtual ~RuntimeComponent() = default;
  /** @brief Runtime borrowed by this component. */
  const ModelRuntime& owning_runtime() const { return runtime_; }
  /** @brief Prepare stable component views before runtime registration is frozen. */
  virtual int seal() = 0;
protected:
  ModelRuntime& runtime_;
};
} // namespace ELFF::Models
