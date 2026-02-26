#pragma once

#include "config/config.hpp"
#include "CommHandle.hpp"
#include "Policy.hpp"

#include <cstddef>
#include <memory>
#include <stdexcept>

#ifdef ELFF_USE_MPI
#include <mpi.h>
#endif 

namespace ELFF {
namespace MPI {

/**
 * @brief Configuration for constructing/initializing a \ref Runtime.
 *
 * This struct controls process-level ownership (MPI init/finalize), behavior
 * when MPI is absent, and how the Runtime stores its base communicator.
 *
 * Typical embedded-library settings (e.g. called from Basilisk-MPI):
 * - mpi_policy = Borrow
 * - no_mpi_policy = Error or SerialFallback
 * - base_comm_policy = DuplicateBase (recommended if library may use
 *   point-to-point communication or wildcard receives)
 */
struct RuntimeConfig
{
  InitPolicy init_policy = InitPolicy::Borrow;
  FallbackPolicy fallback_policy = FallbackPolicy::SerialFallback;
  BaseCommPolicy base_comm_policy = BaseCommPolicy::DuplicateBase;

  /**
   * @brief Whether Runtime::finalize should enforce that no owned child
   * communicator handles remain alive.
   *
   * If true, finalization may fail (e.g. throw/assert) when owned handles
   * created from this Runtime are still active.
   */
  bool require_no_active_handles_on_finalize = true;
};

/**
 * @brief Concrete process-level MPI runtime wrapper for ELFF.
 *
 * @details
 * Runtime models the *process-level* MPI state relevant to the library:
 * whether MPI is available/initialized, who owns MPI initialization, and the
 * library's base communicator (borrowed or duplicated from a caller-provided
 * communicator or MPI_COMM_WORLD).
 *
 * # Key distinction
 * Runtime does **not** represent a "new MPI runtime" when communicators are
 * duplicated. Duplicated communicators are separate objects managed by
 * CommHandle, but they still live under the same process-level MPI
 * initialization.
 *
 * # Intended usage
 * - Usually one Runtime per process (or one per library instance)
 * - Many CommHandle objects created from that Runtime
 * - Embedded mode (Basilisk): Runtime typically borrows MPI initialization
 * - Standalone mode: Runtime may own MPI initialization/finalization
 *
 * # Contracts
 * - If configured with MPIInitPolicy::Borrow, Runtime::finalize must not call
 *   MPI_Finalize.
 * - If configured with BaseCommPolicy::DuplicateBase, Runtime may own and free
 *   a duplicated base communicator.
 * - Owned child communicators (CommHandle) should be destroyed before
 *   Runtime::finalize and always before MPI_Finalize.
 * - Runtime initialization/finalization is explicit; avoid relying on global
 *   static destruction for MPI teardown ordering.
 */
class Runtime
{
public:
  Runtime() noexcept;
  ~Runtime();

  Runtime(const Runtime&) = delete;
  Runtime& operator=(const Runtime&) = delete;
  Runtime(Runtime&&) = delete;
  Runtime& operator=(Runtime&&) = delete;

  /**
   * @brief Initialize using the current configuration and a default base
   * communicator (typically MPI_COMM_WORLD when MPI is enabled).
   *
   * @note This is an explicit lifecycle operation; implementations may make it
   * idempotent, but callers should generally treat repeated init() calls
   * without finalize() as a logic error unless documented otherwise.
   */
  void init();

#ifdef ELFF_USE_MPI
  /**
   * @brief Initialize using a caller-provided base communicator.
   *
   * @param base_comm The caller communicator to borrow or duplicate according
   * to the configured BaseCommPolicy.
   *
   * @warning If duplication is selected, communicator duplication is collective
   *          over \p base_comm.
   */
  void init(MPI_Comm base_comm);
#endif

  /**
   * @brief Finalize the Runtime and release resources it owns.
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
  void finalize();

  /**
   * @brief Set/get runtime configuration.
   *
   * Configuration should generally be set before init().
   */
  void set_config(const RuntimeConfig& cfg);
  const RuntimeConfig& config() const noexcept;

  //! Whether this Runtime has been initialized.
  bool is_initialized() const noexcept;

  //! Whether MPI is enabled/usable for this Runtime instance.
  bool mpi_enabled() const noexcept;

  //! Whether this Runtime owns process-level MPI initialization/finalization.
  bool owns_mpi_runtime() const noexcept;

  //! Whether the base communicator is owned (i.e., duplicated by Runtime).
  bool owns_base_communicator() const noexcept;

  //! Rank in the Runtime's base communicator, or 0 in serial fallback.
  int rank() const noexcept;

  //! Size of the Runtime's base communicator, or 1 in serial fallback.
  int size() const noexcept;

#ifdef ELFF_USE_MPI
  /**
   * @brief Access the Runtime's base communicator.
   *
   * @return MPI_COMM_NULL when MPI is disabled for this Runtime.
   *
   * @note This returns the Runtime base communicator. If the Runtime duplicated
   * the base, this is the duplicated communicator.
   */
  MPI_Comm base_comm() const noexcept;
#endif


#ifdef ELFF_USE_MPI
  /**
   * @brief Create a borrowed handle to the Runtime's base communicator.
   *
   * @return Non-owning communicator handle.
   */
  CommHandle borrow_base_comm() const;
#endif 

#ifdef ELFF_USE_MPI
  /**
   * @brief Create an owned duplicate of the Runtime's base communicator.
   *
   * @return Owning communicator handle with isolated message-matching context.
   *
   * @warning Collective over the Runtime's base communicator.
   */
  CommHandle duplicate_base_comm() const;
#endif 

#ifdef ELFF_USE_MPI
 /**
   * @brief Number of active owned communicator handles tracked by this Runtime.
   *
   * This is primarily a debugging/safety aid used to enforce finalization
   * ordering contracts.
   */
  std::size_t active_owned_handles() const noexcept;
#endif 

private:
  // Runtime may optionally track CommHandle ownership for diagnostics.
  friend class CommHandle;

  //! Internal hooks for owned-handle tracking (used by CommHandle impl).
  void register_owned_handle() noexcept;
  void unregister_owned_handle() noexcept;

  struct Impl;
  std::unique_ptr<Impl> impl_;
};

} // namespace MPI
} // namespace ELFF