#include "error.hpp"
#include "globals.hpp"
// #include "array.hpp"
#include <cstdlib>
#include <iostream>


namespace ELFF
{

#ifdef ELFF_USE_EXCEPTIONS
const char* ErrorException::what() const throw()
{
   return msg.c_str();
}

static ErrorAction elff_error_action = ELFF_ERROR_THROW;
#else
static ErrorAction elff_error_action = ELFF_ERROR_ABORT;
#endif

void set_error_action(ErrorAction action)
{
   // Check if 'action' is valid.
   switch (action)
   {
      case ELFF_ERROR_ABORT: break;
      case ELFF_ERROR_THROW:
#ifdef ELFF_USE_EXCEPTIONS
         break;
#else
         elff_error("set_error_action: ELFF_ERROR_THROW requires the build "
                    "option ELFF_USE_EXCEPTIONS=YES");
         return;
#endif
      default:
         ELFF::err << "\n\nset_error_action: invalid action: " << action
                   << '\n';
         elff_error();
         return;
   }
   elff_error_action = action;
}

ErrorAction get_error_action()
{
   return elff_error_action;
}

namespace internal
{
// defined in globals.cpp
extern bool elff_out_initialized, elff_err_initialized;
}

void elff_backtrace(int mode, int depth)
{

}

void elff_error(const char *msg)
{
   std::ostream &merr = internal::elff_err_initialized ? ELFF::err : std::cerr;
   if (msg)
   {
      // NOTE: By default, each call of the "operator <<" method of the
      // ELFF::err object results in flushing the I/O stream, which can be a
      // very bad thing if all your processors try to do it at the same time.
      merr << "\n\n" << msg << "\n";
   }


#ifdef ELFF_USE_EXCEPTIONS
   if (elff_error_action == ELFF_ERROR_THROW)
   {
      throw ErrorException(msg);
   }
#endif

   std::abort(); // force crash by calling abort
}

void elff_warning(const char *msg)
{
   std::ostream &mout = internal::elff_out_initialized ? ELFF::out : std::cout;
   if (msg)
   {
      mout << "\n\n" << msg << std::endl;
   }
}

}
