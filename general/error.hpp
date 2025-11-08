#pragma once

// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// This file is part of the BEAM library. For more information and source code
// availability visit https://beam.org.
//
// BEAM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

#include "../config/config.hpp"
#include <iomanip>
#include <sstream>

namespace beam
{

/// Action to take when BEAM encounters an error.
enum ErrorAction
{
   BEAM_ERROR_ABORT = 0, /**<
      Abort execution using abort() or MPI_Abort(). This is the default error
      action when the build option BEAM_USE_EXCEPTIONS is set to NO. */
   BEAM_ERROR_THROW      /**<
      Throw an ErrorException. Requires the build option BEAM_USE_EXCEPTIONS=YES
      in which case it is also the default error action. */
};

/// Set the action BEAM takes when an error is encountered.
void set_error_action(ErrorAction action);
/// Get the action BEAM takes when an error is encountered.
ErrorAction get_error_action();

#ifdef BEAM_USE_EXCEPTIONS
/** @brief Exception class thrown when BEAM encounters an error and the current
    ErrorAction is set to BEAM_ERROR_THROW. */
class ErrorException: public std::exception
{
private:
   std::string msg;
public:
   explicit ErrorException(const std::string & in_msg) : msg(in_msg) { }
   virtual ~ErrorException() throw() { }
   virtual const char* what() const throw();
};
#endif

void beam_backtrace(int mode = 0, int depth = -1);

/** @brief Function called when an error is encountered. Used by the macros
    BEAM_ABORT, BEAM_ASSERT, BEAM_VERIFY. */
[[noreturn]]
void beam_error(const char *msg = NULL);

/// Function called by the macro BEAM_WARNING.

void beam_warning(const char *msg = NULL);
}

#ifndef _BEAM_FUNC_NAME
#ifndef _MSC_VER
// This is nice because it shows the class and method name
#define _BEAM_FUNC_NAME __PRETTY_FUNCTION__
// This one is C99 standard.
//#define _BEAM_FUNC_NAME __func__
#else
// for Visual Studio C++
#define _BEAM_FUNC_NAME __FUNCSIG__
#endif
#endif

#define BEAM_LOCATION \
   "\n ... in function: " << _BEAM_FUNC_NAME << \
   "\n ... in file: " << __FILE__ << ':' << __LINE__ << '\n'

// Common error message and abort macro
#define _BEAM_MESSAGE(msg, fn)                                          \
   {                                                                    \
      std::ostringstream beamMsgStream;                                 \
      beamMsgStream << std::setprecision(16);                           \
      beamMsgStream << std::setiosflags(std::ios_base::scientific);     \
      beamMsgStream << msg << BEAM_LOCATION;                            \
      beam::fn(beamMsgStream.str().c_str());                            \
   }

// Outputs lots of useful information and aborts.
// For all of these functions, "msg" is pushed to an ostream, so you can
// write useful (if complicated) error messages instead of writing
// out to the screen first, then calling abort.  For example:
// BEAM_ABORT( "Unknown geometry type: " << type );
#define BEAM_ABORT(msg) _BEAM_MESSAGE("BEAM abort: " << msg, beam_error)

// Does a check, and then outputs lots of useful information if the test fails
#define BEAM_VERIFY(x, msg)                             \
   if (!(x))                                            \
   {                                                    \
      _BEAM_MESSAGE("Verification failed: ("            \
                    << #x << ") is false:\n --> " << msg, beam_error); \
   }

// Use this if the only place your variable is used is in ASSERTs
// For example, this code snippet:
//   int err = MPI_Reduce(ldata, maxdata, 5, MPI_INT, MPI_MAX, 0, MyComm);
//   BEAM_CONTRACT_VAR(err);
//   BEAM_ASSERT( err == 0, "MPI_Reduce gave an error with length "
//                       << ldata );
#define BEAM_CONTRACT_VAR(x) (void)(x)

// Now set up some optional checks, but only if the right flags are on
#ifdef BEAM_DEBUG

#define BEAM_ASSERT(x, msg)                             \
   if (!(x))                                            \
   {                                                    \
      _BEAM_MESSAGE("Assertion failed: ("               \
                    << #x << ") is false:\n --> " << msg, beam_error); \
   }

// A macro that exposes its argument in debug mode only.
#define BEAM_DEBUG_DO(x) x

#else

// Get rid of all this code, since we're not checking.
#define BEAM_ASSERT(x, msg)

// A macro that exposes its argument in debug mode only.
#define BEAM_DEBUG_DO(x)

#endif

// Generate a warning message - always generated, regardless of BEAM_DEBUG.
#define BEAM_WARNING(msg) _BEAM_MESSAGE("BEAM Warning: " << msg, beam_warning)

// Macro that checks (in BEAM_DEBUG mode) that i is in the range [imin,imax).
#define BEAM_ASSERT_INDEX_IN_RANGE(i,imin,imax) \
   BEAM_ASSERT((imin) <= (i) && (i) < (imax), \
   "invalid index " #i << " = " << (i) << \
   ", valid range is [" << (imin) << ',' << (imax) << ')')

