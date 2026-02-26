#include "Runtime.hpp"

#include <stdexcept> // std::logic_error, std::runtime_error
#include <utility>   // std::move
#include <string>    // std::string, std::to_string

namespace ELFF {
namespace MPI {

namespace {

#ifdef ELFF_USE_MPI
/**
 * @brief Throw a formatted C++ exception for an MPI API failure.
 *
 * @param where Human-readable description of the failing MPI call/site.
 * @param errcode MPI error code returned by the MPI routine.
 *
 * @return This function does not return; it always throws.
 *
 * @throws std::runtime_error Always thrown with a formatted message.
 */
inline void
throw_mpi_error(const char* where, int errcode)
{
  throw std::runtime_error(std::string(where) + " failed with MPI error code " +
                           std::to_string(errcode));
}
#endif

#ifdef ELFF_USE_MPI
/**
 * @brief Query whether the MPI runtime has been initialized in this process.
 *
 * @return true if MPI_Init/MPI_Init_thread has been called, false otherwise.
 */
inline bool
mpi_is_initialized()
{
  int flag = 0;
  MPI_Initialized(&flag);
  return flag != 0;
}
#endif

#ifdef ELFF_USE_MPI
/**
 * @brief Query whether the MPI runtime has already been finalized.
 *
 * @return true if MPI_Finalize has been called, false otherwise.
 */
inline bool
mpi_is_finalized()
{
  int flag = 0;
  MPI_Finalized(&flag);
  return flag != 0;
}
#endif

} // namespace

/**
 * @brief Private implementation state for @ref Runtime (PIMPL idiom).
 *
 * @details
 * Stores configuration, lifecycle state, rank/size information, and (in MPI
 * builds) the base communicator handle plus optional diagnostics for tracking
 * owned child communicator handles.
 */
struct Runtime::Impl
{
  RuntimeConfig config{};

  bool initialized = false;
  bool mpi_enabled = false;
  bool owns_mpi_runtime = false;

  int rank = 0;
  int size = 1;

#ifdef ELFF_USE_MPI
  CommHandle base_comm_handle{};
  std::size_t active_owned_handles =
    0; // NOTE: requires CommHandle callbacks to be fully accurate
#endif
};

/**
 * @brief Construct a Runtime in the uninitialized state.
 *
 * @details
 * Allocates the internal implementation object and initializes the runtime as
 * disabled/uninitialized until @ref Runtime::init is called.
 */
Runtime::Runtime() noexcept
  : impl_(std::make_unique<Impl>())
{
}

/**
 * @brief Destroy the Runtime and perform best-effort cleanup.
 *
 * @details
 * If the runtime is still initialized, this destructor attempts to call
 * @ref Runtime::finalize with strict active-handle enforcement temporarily
 * disabled. Exceptions are swallowed to preserve noexcept destructor behavior.
 *
 * @note Prefer explicit @ref Runtime::finalize calls in user code for clearer
 *       teardown ordering and error handling.
 */
Runtime::~Runtime()
{
  if (!impl_)
    return;

  // Best-effort cleanup without throwing from destructor.
  if (impl_->initialized) {
    const bool saved = impl_->config.require_no_active_handles_on_finalize;
    impl_->config.require_no_active_handles_on_finalize = false;
    try {
      finalize();
    } catch (...) {
      // swallow in destructor
    }
    impl_->config.require_no_active_handles_on_finalize = saved;
  }
}

#ifdef ELFF_USE_MPI
/**
 * @brief (MPI) Initialize using the current configuration and a default base
 * communicator (typically MPI_COMM_WORLD when MPI is enabled).
 *
 * @note This is an explicit lifecycle operation; implementations may make it
 * idempotent, but callers should generally treat repeated init() calls
 * without finalize() as a logic error unless documented otherwise.
 */
void
Runtime::init()
{
  init(
    MPI_COMM_NULL); // MPI_COMM_NULL => use default (WORLD) once MPI is active
}
#else
/**
 * @brief (SERIAL) Initialize the runtime using the current configuration.
 *
 * @details
 * In a non-MPI build, initialization succeeds only if serial fallback is
 * permitted by the configured fallback policy. The runtime is then initialized
 * in serial mode with rank=0 and size=1.
 *
 * @note This is an explicit lifecycle operation; callers should generally treat
 * repeated init() calls without finalize() as a logic error.
 *
 * @throws std::logic_error If the runtime is already initialized.
 * @throws std::runtime_error If serial fallback is disallowed.
 */
void
Runtime::init()
{
  if (impl_->initialized) {
    throw std::logic_error(
      "ELFF::MPI::Runtime::init: runtime already initialized");
  }

  if (impl_->config.fallback_policy == FallbackPolicy::Error) {
    throw std::runtime_error("ELFF::MPI::Runtime::init: ELFF built without MPI "
                             "and fallback policy is Error");
  }

  // Serial fallback
  impl_->initialized = true;
  impl_->mpi_enabled = false;
  impl_->owns_mpi_runtime = false;
  impl_->rank = 0;
  impl_->size = 1;
}
#endif

#ifdef ELFF_USE_MPI
/**
 * @brief (MPI) Initialize using a caller-provided base communicator.
 *
 * @param base_comm The caller communicator to borrow or duplicate according
 * to the configured BaseCommPolicy.
 *
 * @warning If duplication is selected, communicator duplication is collective
 *          over \p base_comm.
 */
void
Runtime::init(MPI_Comm base_comm)
{
  if (impl_->initialized) {
    throw std::logic_error(
      "ELFF::MPI::Runtime::init: runtime already initialized");
  }

  // Reset to a known state before attempting initialization.
  impl_->initialized = false;
  impl_->mpi_enabled = false;
  impl_->owns_mpi_runtime = false;
  impl_->rank = 0;
  impl_->size = 1;
  impl_->base_comm_handle = CommHandle{};
  impl_->active_owned_handles = 0;

  bool initialized = mpi_is_initialized();
  if (!initialized) {
    if (impl_->config.init_policy == InitPolicy::Own) {
      const int err = MPI_Init(nullptr, nullptr);
      if (err != MPI_SUCCESS) {
        throw_mpi_error("ELFF::MPI::Runtime::init (MPI_Init)", err);
      }
      impl_->owns_mpi_runtime = true;
      initialized = true;
    } else {
      // Borrow policy, but MPI not initialized.
      if (impl_->config.fallback_policy == FallbackPolicy::Error) {
        throw std::runtime_error(
          "ELFF::MPI::Runtime::init: MPI is not initialized in Borrow mode");
      }

      // Serial fallback.
      impl_->initialized = true;
      impl_->mpi_enabled = false;
      impl_->owns_mpi_runtime = false;
      impl_->rank = 0;
      impl_->size = 1;
      return;
    }
  }

  if (mpi_is_finalized()) {
    throw std::runtime_error(
      "ELFF::MPI::Runtime::init: MPI has already been finalized");
  }

  // Resolve default communicator now that MPI is active.
  MPI_Comm source = (base_comm != MPI_COMM_NULL) ? base_comm : MPI_COMM_WORLD;

  if (impl_->config.base_comm_policy == BaseCommPolicy::DuplicateBase) {
    // Collective over source.
    impl_->base_comm_handle = CommHandle::duplicate(source);
  } else {
    impl_->base_comm_handle = CommHandle::borrow(source);
  }

  if (!impl_->base_comm_handle) {
    throw std::runtime_error("ELFF::MPI::Runtime::init: failed to establish a "
                             "valid base communicator");
  }

  const MPI_Comm comm = impl_->base_comm_handle.get();
  {
    const int err_rank = MPI_Comm_rank(comm, &impl_->rank);
    if (err_rank != MPI_SUCCESS) {
      throw_mpi_error("ELFF::MPI::Runtime::init (MPI_Comm_rank)", err_rank);
    }
    const int err_size = MPI_Comm_size(comm, &impl_->size);
    if (err_size != MPI_SUCCESS) {
      throw_mpi_error("ELFF::MPI::Runtime::init (MPI_Comm_size)", err_size);
    }
  }

  impl_->mpi_enabled = true;
  impl_->initialized = true;
}
#endif

#ifdef ELFF_USE_MPI
/**
 * @brief (MPI) Finalize the Runtime and release resources it owns.
 *
 * This may include:
 * - freeing an owned duplicated base communicator
 * - (only in MPIInitPolicy::Own mode) calling MPI_Finalize
 *
 * This must not free communicators borrowed from the caller.
 *
 * @warning If active owned CommHandle objects still exist, finalization
 * may fail or assert depending on implementation and configuration.
 */
void
Runtime::finalize()
{
  if (!impl_->initialized) {
    return;
  }

  if (impl_->config.require_no_active_handles_on_finalize &&
      impl_->active_owned_handles != 0) {
    throw std::logic_error(
      "ELFF::MPI::Runtime::finalize: active owned CommHandle instances remain");
  }

  // Release base communicator first (before possibly finalizing MPI).
  // Use move-assignment to a default handle so owned communicators are freed
  // using CommHandle's destructor/move-assignment semantics.
  impl_->base_comm_handle = CommHandle{};

  if (impl_->owns_mpi_runtime) {
    if (mpi_is_initialized() && !mpi_is_finalized()) {
      const int err = MPI_Finalize();
      if (err != MPI_SUCCESS) {
        throw_mpi_error("ELFF::MPI::Runtime::finalize (MPI_Finalize)", err);
      }
    }
    impl_->owns_mpi_runtime = false;
  }

  impl_->mpi_enabled = false;
  impl_->rank = 0;
  impl_->size = 1;
  impl_->initialized = false;
}
#else
/**
 * @brief (SERIAL) Finalize the Runtime and reset it to the uninitialized state.
 *
 * @details
 * In non-MPI builds, this clears serial-runtime state only. No MPI resources
 * are released and no MPI finalization is attempted.
 */
void
Runtime::finalize()
{
  if (!impl_->initialized) {
    return;
  }

  impl_->mpi_enabled = false;
  impl_->rank = 0;
  impl_->size = 1;
  impl_->initialized = false;
}

#endif

/**
 * @brief Set/get runtime configuration.
 *
 * Configuration should generally be set before init().
 */
void
Runtime::set_config(const RuntimeConfig& cfg)
{
  if (impl_->initialized) {
    throw std::logic_error(
      "ELFF::MPI::Runtime::set_config: cannot change config after init()");
  }
  impl_->config = cfg;
}

/**
 * @brief Get the current runtime configuration.
 *
 * @return Const reference to the stored @ref RuntimeConfig.
 */
const RuntimeConfig&
Runtime::config() const noexcept
{
  return impl_->config;
}

/**
 * @brief Check whether this Runtime has been initialized.
 *
 * @return true if @ref init has completed successfully and @ref finalize has
 *         not yet been called; false otherwise.
 */
bool
Runtime::is_initialized() const noexcept
{
  return impl_->initialized;
}

/**
 * @brief Check whether MPI is enabled/usable for this Runtime instance.
 *
 * @details
 * This may be false in non-MPI builds or when the runtime is operating in
 * serial fallback mode.
 *
 * @return true if MPI is active for this Runtime; false otherwise.
 */
bool
Runtime::mpi_enabled() const noexcept
{
  return impl_->mpi_enabled;
}

/**
 * @brief Check whether this Runtime owns process-level MPI initialization.
 *
 * @details
 * When true, the runtime is responsible for calling MPI_Finalize during
 * @ref finalize. When false, MPI lifetime is managed by the caller/host.
 *
 * @return true if this Runtime owns MPI init/finalize responsibility.
 */
bool
Runtime::owns_mpi_runtime() const noexcept
{
  return impl_->owns_mpi_runtime;
}

#ifdef ELFF_USE_MPI
/**
 * @brief (MPI) Check whether the Runtime owns its base communicator.
 *
 * @details
 * This is true when the base communicator was duplicated by the Runtime (e.g.
 * BaseCommPolicy::DuplicateBase), and false when the Runtime is borrowing a
 * caller-provided communicator.
 *
 * @return true if the base communicator is owned by the Runtime.
 */
bool
Runtime::owns_base_communicator() const noexcept
{
  return impl_->base_comm_handle.owns_communicator();
}
#else
/**
 * @brief (SERIAL) Report base communicator ownership in a non-MPI build.
 *
 * @return Always false, since no MPI communicator exists in serial builds.
 */
bool
Runtime::owns_base_communicator() const noexcept
{
  return false;
}
#endif

/**
 * @brief Rank in the Runtime's base communicator, or 0 in serial fallback.
 */
int
Runtime::rank() const noexcept
{
  return impl_->rank;
}

/**
 * @brief Size of the Runtime's base communicator, or 1 in serial fallback.
 */
int
Runtime::size() const noexcept
{
  return impl_->size;
}

#ifdef ELFF_USE_MPI
/**
 * @brief (MPI) Access the Runtime's base communicator.
 *
 * @return MPI_COMM_NULL when MPI is disabled for this Runtime.
 *
 * @note This returns the Runtime base communicator. If the Runtime duplicated
 * the base, this is the duplicated communicator.
 */
MPI_Comm
Runtime::base_comm() const noexcept
{
  if (!impl_->initialized || !impl_->mpi_enabled)
    return MPI_COMM_NULL;
  return impl_->base_comm_handle.get();
}
#endif

#ifdef ELFF_USE_MPI
/**
 * @brief (MPI) Create a borrowed handle to the Runtime's base communicator.
 *
 * @return Non-owning communicator handle.
 */
CommHandle
Runtime::borrow_base_comm() const
{
  if (!impl_->initialized) {
    throw std::logic_error(
      "ELFF::MPI::Runtime::borrow_base_comm: runtime is not initialized");
  }

  if (!impl_->mpi_enabled) {
    return CommHandle{}; // serial fallback
  }

  return CommHandle::borrow(impl_->base_comm_handle.get());
}
#endif

#ifdef ELFF_USE_MPI
/**
 * @brief (MPI) Create an owned duplicate of the Runtime's base communicator.
 *
 * @return Owning communicator handle with isolated message-matching context.
 *
 * @warning Collective over the Runtime's base communicator.
 */
CommHandle
Runtime::duplicate_base_comm() const
{
  if (!impl_->initialized) {
    throw std::logic_error(
      "ELFF::MPI::Runtime::duplicate_base_comm: runtime is not initialized");
  }

  if (!impl_->mpi_enabled) {
    return CommHandle{}; // serial fallback
  }

  // NOTE: active_owned_handles tracking is only accurate if CommHandle
  // implementation is extended to call register/unregister hooks.
  return CommHandle::duplicate(impl_->base_comm_handle.get());
}
#endif

#ifdef ELFF_USE_MPI
/**
 * @brief (MPI) Number of active owned communicator handles tracked by this Runtime.
 *
 * This is primarily a debugging/safety aid used to enforce finalization
 * ordering contracts.
 */
std::size_t
Runtime::active_owned_handles() const noexcept
{
  return impl_->active_owned_handles;
}
#endif

#ifdef ELFF_USE_MPI
/**
 * @brief (MPI) Increment the internal count of owned child communicator handles.
 *
 * @details
 * This is an internal diagnostic/safety hook intended to be called by
 * @ref CommHandle when a Runtime-tracked owned communicator is created.
 *
 * @return None.
 */
void
Runtime::register_owned_handle() noexcept
{
  ++impl_->active_owned_handles;
}
#else 
/**
 * @brief (SERIAL) No-op owned-handle registration hook.
 *
 * @details
 * In non-MPI builds there are no MPI communicator handles to track, so this
 * function intentionally does nothing.
 *
 * @return None.
 */
void
Runtime::register_owned_handle() noexcept {}
#endif

#ifdef ELFF_USE_MPI
/**
 * @brief (MPI) Decrement the internal count of owned child communicator handles.
 *
 * @details
 * This is an internal diagnostic/safety hook intended to be called by
 * @ref CommHandle when a Runtime-tracked owned communicator is destroyed or
 * released. The counter is clamped at zero to avoid underflow.
 *
 * @return None.
 */
void
Runtime::unregister_owned_handle() noexcept
{
  if (impl_->active_owned_handles > 0)
    --impl_->active_owned_handles;
}
#else
/**
 * @brief (SERIAL) No-op owned-handle unregistration hook.
 *
 * @details
 * In non-MPI builds there are no MPI communicator handles to track, so this
 * function intentionally does nothing.
 *
 * @return None.
 */
void
Runtime::unregister_owned_handle() noexcept
{}
#endif

} // namespace MPI
} // namespace ELFF