#include "CommHandle.hpp"

#include <stdexcept> // std::invalid_argument, std::runtime_error
#include <string>    // std::string, std::to_string
#include <utility>   // std::move

namespace ELFF {
namespace MPI {

namespace {

#ifdef ELFF_USE_MPI
/**
 * @brief Throw a runtime error describing an MPI API failure.
 *
 * @param where Human-readable description of the call site / MPI operation.
 * @param errcode MPI error code returned by the failed MPI routine.
 *
 * @return This function does not return; it always throws.
 *
 * @throws std::runtime_error Always thrown with a formatted error message.
 */
inline void
throw_mpi_error(const char* where, int errcode)
{
  // Keep this simple; MPI's default error handler often aborts anyway.
  throw std::runtime_error(std::string(where) + " failed with MPI error code " +
                           std::to_string(errcode));
}
#endif

} // namespace

/**
 * @brief Construct an invalid/null communicator handle.
 *
 * @details
 * The default-constructed handle is non-owning and evaluates to false.
 *
 * @return A default-initialized invalid handle (constructor).
 */
#ifdef ELFF_USE_MPI
CommHandle::CommHandle() noexcept
  : comm_(MPI_COMM_NULL)
  , owned_(false)
  , valid_(false)
{
}
#else
CommHandle::CommHandle() noexcept
  : owned_(false)
  , valid_(false)
{
}
#endif

/**
 * @brief Construct a handle from explicit ownership and communicator state.
 *
 * @details
 * This constructor is intended for internal use by factory methods such as
 * borrow() and duplicate().
 *
 * @param owned True if the handle owns the communicator and should free it.
 * @param comm The communicator to wrap (MPI build only).
 *
 * @return A handle initialized from the given ownership/communicator state
 *         (constructor).
 */
#ifdef ELFF_USE_MPI
CommHandle::CommHandle(bool owned, MPI_Comm comm) noexcept
  : comm_(comm)
  , owned_(owned)
  , valid_(comm != MPI_COMM_NULL)
{
}
#else
CommHandle::CommHandle(bool owned) noexcept
  : owned_(owned)
  , valid_(false)
{
}
#endif

/**
 * @brief Destroy the communicator handle.
 *
 * @details
 * If this handle owns a valid communicator and MPI is initialized (and not yet
 * finalized), the communicator is freed via MPI_Comm_free(). If MPI is already
 * finalized (or was never initialized), the handle is simply invalidated
 * locally.
 *
 * @return None.
 */
#ifdef ELFF_USE_MPI
CommHandle::~CommHandle()
{
  // Destructor semantics: free only if this handle owns a valid communicator.
  if (owned_ && valid_ && comm_ != MPI_COMM_NULL) {
    // Be defensive: MPI_Comm_free is only valid while MPI is initialized
    // and not yet finalized.
    int is_initialized = 0;
    int is_finalized = 0;
    MPI_Initialized(&is_initialized);
    if (is_initialized) {
      MPI_Finalized(&is_finalized);
      if (!is_finalized) {
        // MPI_Comm_free nulls the communicator on success.
        (void)MPI_Comm_free(&comm_);
      } else {
        // MPI already finalized; cannot safely free. Invalidate locally.
        comm_ = MPI_COMM_NULL;
      }
    } else {
      // MPI was never initialized; cannot safely free. Invalidate locally.
      comm_ = MPI_COMM_NULL;
    }
  }
  owned_ = false;
  valid_ = false;
}
#else
CommHandle::~CommHandle()
{
  owned_ = false;
  valid_ = false;
}
#endif

/**
 * @brief Move-construct a communicator handle.
 *
 * @param other Source handle whose state is transferred into this object.
 *
 * @return A new handle that takes ownership/validity state from @p other
 *         (constructor).
 */
#ifdef ELFF_USE_MPI
CommHandle::CommHandle(CommHandle&& other) noexcept
  : comm_(other.comm_)
  , owned_(other.owned_)
  , valid_(other.valid_)
{
  other.comm_ = MPI_COMM_NULL;
  other.owned_ = false;
  other.valid_ = false;
}
#else
CommHandle::CommHandle(CommHandle&& other) noexcept
  : owned_(other.owned_)
  , valid_(other.valid_)
{
  other.owned_ = false;
  other.valid_ = false;
}
#endif

/**
 * @brief Move-assign from another communicator handle.
 *
 * @details
 * Any currently owned communicator held by this object is released first
 * (subject to MPI being initialized and not finalized), then the state from
 * @p other is transferred and @p other is invalidated.
 *
 * @param other Source handle whose state is transferred into this object.
 *
 * @return CommHandle& Reference to `*this`.
 */
#ifdef ELFF_USE_MPI
CommHandle&
CommHandle::operator=(CommHandle&& other) noexcept
{
  if (this == &other)
    return *this;

  // First, free our currently-owned communicator (destructor semantics).
  if (owned_ && valid_ && comm_ != MPI_COMM_NULL) {
    int is_initialized = 0;
    int is_finalized = 0;
    MPI_Initialized(&is_initialized);
    if (is_initialized) {
      MPI_Finalized(&is_finalized);
      if (!is_finalized) {
        (void)MPI_Comm_free(&comm_);
      } else {
        comm_ = MPI_COMM_NULL;
      }
    } else {
      comm_ = MPI_COMM_NULL;
    }
  }

  // Steal state
  comm_ = other.comm_;
  owned_ = other.owned_;
  valid_ = other.valid_;

  other.comm_ = MPI_COMM_NULL;
  other.owned_ = false;
  other.valid_ = false;

  return *this;
}
#else
CommHandle&
CommHandle::operator=(CommHandle&& other) noexcept
{
  if (this == &other)
    return *this;

  // No MPI build: just move flags.
  owned_ = other.owned_;
  valid_ = other.valid_;

  other.owned_ = false;
  other.valid_ = false;

  return *this;
}
#endif

#ifdef ELFF_USE_MPI
/**
 * @brief Create a non-owning communicator handle.
 *
 * @param comm Communicator to wrap. May be MPI_COMM_NULL, in which case the
 *             returned handle is invalid/non-owning.
 *
 * @return CommHandle A borrowed (non-owning) communicator handle.
 */
CommHandle
CommHandle::borrow(MPI_Comm comm)
{
  // Borrowing MPI_COMM_NULL is allowed and yields an invalid handle.
  return CommHandle(false, comm);
}
#endif

#ifdef ELFF_USE_MPI
/**
 * @brief Duplicate a communicator and return an owning handle.
 *
 * @details
 * This performs MPI_Comm_dup(), which is collective over the source
 * communicator. The returned handle owns the duplicated communicator and will
 * free it in its destructor (subject to MPI lifetime constraints).
 *
 * @param comm Source communicator to duplicate. Must not be MPI_COMM_NULL.
 *
 * @return CommHandle An owning handle wrapping the duplicated communicator.
 *
 * @throws std::invalid_argument If @p comm is MPI_COMM_NULL.
 * @throws std::runtime_error If MPI is not initialized, already finalized, or
 *         if MPI_Comm_dup() fails.
 */
CommHandle
CommHandle::duplicate(MPI_Comm comm)
{
  if (comm == MPI_COMM_NULL) {
    throw std::invalid_argument(
      "ELFF::MPI::CommHandle::duplicate: comm is MPI_COMM_NULL");
  }

  int is_initialized = 0;
  int is_finalized = 0;

  MPI_Initialized(&is_initialized);
  if (!is_initialized) {
    throw std::runtime_error(
      "ELFF::MPI::CommHandle::duplicate: MPI is not initialized");
  }

  MPI_Finalized(&is_finalized);
  if (is_finalized) {
    throw std::runtime_error(
      "ELFF::MPI::CommHandle::duplicate: MPI is already finalized");
  }

  MPI_Comm dup = MPI_COMM_NULL;
  const int err = MPI_Comm_dup(comm, &dup);
  if (err != MPI_SUCCESS) {
    throw_mpi_error("ELFF::MPI::CommHandle::duplicate (MPI_Comm_dup)", err);
  }

  return CommHandle(true, dup);
}
#endif

#ifdef ELFF_USE_MPI
/**
 * @brief Get the underlying MPI communicator.
 *
 * @return MPI_Comm The wrapped communicator, possibly MPI_COMM_NULL.
 */
MPI_Comm
CommHandle::get() const noexcept
{
  return comm_;
}
#endif

/**
 * @brief Test whether the handle currently wraps a valid communicator.
 *
 * @return bool True if the handle is valid; false otherwise.
 */
CommHandle::operator bool() const noexcept
{
  return valid_;
}

/**
 * @brief Query whether this handle owns the communicator.
 *
 * @return bool True if this handle will free the communicator; false if it is
 *         borrowed or invalid.
 */
bool
CommHandle::owns_communicator() const noexcept
{
  return owned_;
}

/**
 * @brief Release ownership (if any) and invalidate the handle without freeing.
 *
 * @details
 * This is intended for advanced ownership-transfer scenarios. After reset(),
 * the handle becomes invalid and non-owning. In MPI builds, the stored
 * communicator is set to MPI_COMM_NULL but not freed.
 *
 * @return None.
 */
#ifdef ELFF_USE_MPI
void
CommHandle::reset() noexcept
{
  comm_ = MPI_COMM_NULL;
  owned_ = false;
  valid_ = false;
}
#else
void
CommHandle::reset() noexcept
{
  owned_ = false;
  valid_ = false;
}
#endif

} // namespace MPI
} // namespace ELFF