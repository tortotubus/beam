#pragma once

namespace ELFF {
namespace MPI {

/**
 * @brief Policy for MPI process-level initialization ownership.
 *
 * This policy applies to a \ref Runtime object (not to communicators).
 *
 * - Borrow: the host application (e.g. Basilisk, another library, or main())
 *   is responsible for calling MPI_Init/MPI_Init_thread and MPI_Finalize.
 * - Own: the Runtime may initialize/finalize MPI itself (typically only for
 *   standalone executables, not embedded-library use).
 *
 * @note MPI initialization/finalization is process-level and happens at most
 *       once per process. This is distinct from communicator duplication/free.
 */
enum class InitPolicy
{
  Borrow,
  Own
};

/**
 * @brief Behavior when MPI is unavailable or not yet initialized.
 *
 * This policy determines how \ref Runtime::init behaves when the library is
 * built with MPI support but the process has not yet called MPI_Init.
 *
 * - Error: fail initialization (e.g. throw) rather than proceeding.
 * - SerialFallback: treat the runtime as disabled/serial (rank=0,size=1).
 */
enum class FallbackPolicy
{
  Error,
  SerialFallback
};

/**
 * @brief Policy for how a Runtime should hold its base communicator.
 *
 * - BorrowBase: store/use the caller communicator directly (non-owning).
 * - DuplicateBase: duplicate the caller communicator to obtain an isolated
 *   message-matching namespace for library traffic.
 *
 * @warning Communicator duplication is collective over the source communicator.
 */
enum class BaseCommPolicy
{
  BorrowBase,
  DuplicateBase
};

} // namespace MPI
} // namespace ELFF