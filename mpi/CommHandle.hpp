#pragma once

#include "config/config.hpp"

#ifdef ELFF_USE_MPI
#include <mpi.h>
#endif

namespace ELFF {
namespace MPI {

/**
 * @brief RAII wrapper for an MPI communicator handle.
 *
 * @details
 * A CommHandle represents a communicator object that is either:
 * - borrowed (non-owning): wraps an existing communicator without freeing it
 * - owned (owning): typically created via MPI_Comm_dup (or split) and frees
 *   itself in the destructor
 *
 * This class is intentionally move-only to make ownership explicit.
 *
 * ## Design intent
 * - CommHandle manages *communicator lifetime* (MPI_Comm_dup/free),
 *   not *MPI runtime lifetime* (MPI_Init/finalize).
 * - It is safe to have many CommHandle instances under one Runtime.
 * - Borrowed handles never free the wrapped communicator.
 *
 * ## Lifetime contract
 * - Any owned communicator must be destroyed (or released) before MPI_Finalize.
 * - Therefore, owned CommHandle objects should not outlive the Runtime
 *   that conceptually governs their use, and must not survive past
 * process-level MPI finalization.
 *
 * ## Collective contract
 * - Factory operations that duplicate/split communicators are collective over
 *   the source communicator and must be called consistently by participating
 *   ranks.
 */
class CommHandle
{
public:
  /**
   * @brief Construct an invalid/null communicator handle.
   *
   * @details
   * The default-constructed handle is non-owning and evaluates to false.
   *
   * @return A default-initialized invalid handle (constructor).
   */
  CommHandle() noexcept;

  /**
   * @brief Destroy the communicator handle.
   *
   * @details
   * If this handle owns a valid communicator and MPI is initialized (and not
   * yet finalized), the communicator is freed via MPI_Comm_free(). If MPI is
   * already finalized (or was never initialized), the handle is simply
   * invalidated locally.
   *
   * @return None.
   */
  ~CommHandle();

  // Move-only ownership semantics.
  CommHandle(CommHandle&& other) noexcept;
  CommHandle& operator=(CommHandle&& other) noexcept;

  CommHandle(const CommHandle&) = delete;
  CommHandle& operator=(const CommHandle&) = delete;

#ifdef ELFF_USE_MPI
  /**
   * @brief Create a non-owning wrapper around an existing communicator.
   *
   * @param comm Communicator to wrap (may be MPI_COMM_NULL for an invalid
   * handle).
   * @return A borrowed handle.
   *
   * @note The caller remains responsible for the communicator lifetime.
   */
  static CommHandle borrow(MPI_Comm comm);
#else
#endif

#ifdef ELFF_USE_MPI
  /**
   * @brief Duplicate a communicator and return an owning handle.
   *
   * @param comm Source communicator to duplicate.
   * @return An owning handle to the duplicated communicator.
   *
   * @warning Collective over @p comm.
   * @warning Intended for use only while MPI is initialized.
   */
  static CommHandle duplicate(MPI_Comm comm);
#else
#endif

#ifdef ELFF_USE_MPI
  /**
   * @brief Access the underlying MPI communicator.
   *
   * @return The wrapped communicator (possibly MPI_COMM_NULL).
   */
  MPI_Comm get() const noexcept;
#else
#endif

  /**
   * @brief Whether this handle currently wraps a valid communicator.
   *
   * In non-MPI builds this will typically be false.
   */
  explicit operator bool() const noexcept;

  /**
   * @brief Whether this handle owns the communicator and will free it.
   */
  bool owns_communicator() const noexcept;

  /**
   * @brief Release ownership without freeing and reset to invalid.
   *
   * Intended only for advanced ownership handoff scenarios.
   */
  void reset() noexcept;

private:
#ifdef ELFF_USE_MPI
  /**
   * @brief Internal constructor used by factory methods.
   *
   * Semantics are implementation-defined; exposed here only to keep the public
   * API minimal and explicit.
   */
  CommHandle(bool owned, MPI_Comm comm) noexcept;
#else
  CommHandle(bool owned) noexcept;
#endif

#ifdef ELFF_USE_MPI
  MPI_Comm comm_ = MPI_COMM_NULL;
#endif
  bool owned_ = false;
  bool valid_ = false;
};

} // namespace MPI
} // namespace ELFF
