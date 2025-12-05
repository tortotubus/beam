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

#include "error.hpp"
#include "globals.hpp"
// #include "array.hpp"
#include <cstdlib>
#include <iostream>


namespace beam
{

#ifdef BEAM_USE_EXCEPTIONS
const char* ErrorException::what() const throw()
{
   return msg.c_str();
}

static ErrorAction beam_error_action = BEAM_ERROR_THROW;
#else
static ErrorAction beam_error_action = BEAM_ERROR_ABORT;
#endif

void set_error_action(ErrorAction action)
{
   // Check if 'action' is valid.
   switch (action)
   {
      case BEAM_ERROR_ABORT: break;
      case BEAM_ERROR_THROW:
#ifdef BEAM_USE_EXCEPTIONS
         break;
#else
         beam_error("set_error_action: BEAM_ERROR_THROW requires the build "
                    "option BEAM_USE_EXCEPTIONS=YES");
         return;
#endif
      default:
         beam::err << "\n\nset_error_action: invalid action: " << action
                   << '\n';
         beam_error();
         return;
   }
   beam_error_action = action;
}

ErrorAction get_error_action()
{
   return beam_error_action;
}

namespace internal
{
// defined in globals.cpp
extern bool beam_out_initialized, beam_err_initialized;
}

void beam_backtrace(int mode, int depth)
{

}

void beam_error(const char *msg)
{
   std::ostream &merr = internal::beam_err_initialized ? beam::err : std::cerr;
   if (msg)
   {
      // NOTE: By default, each call of the "operator <<" method of the
      // beam::err object results in flushing the I/O stream, which can be a
      // very bad thing if all your processors try to do it at the same time.
      merr << "\n\n" << msg << "\n";
   }


#ifdef BEAM_USE_EXCEPTIONS
   if (beam_error_action == BEAM_ERROR_THROW)
   {
      throw ErrorException(msg);
   }
#endif

   std::abort(); // force crash by calling abort
}

void beam_warning(const char *msg)
{
   std::ostream &mout = internal::beam_out_initialized ? beam::out : std::cout;
   if (msg)
   {
      mout << "\n\n" << msg << std::endl;
   }
}

}
