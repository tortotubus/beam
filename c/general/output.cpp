#include "elff/c/general/output.h"

#include "elff/general/globals.hpp"

extern "C" {

void elff_set_out_prefix(const char *prefix)
{
  if (prefix == nullptr) {
    ELFF::ClearOutPrefix();
    return;
  }

  ELFF::SetOutPrefix(prefix);
}

void elff_clear_out_prefix(void)
{
  ELFF::ClearOutPrefix();
}

void elff_set_err_prefix(const char *prefix)
{
  if (prefix == nullptr) {
    ELFF::ClearErrPrefix();
    return;
  }

  ELFF::SetErrPrefix(prefix);
}

void elff_clear_err_prefix(void)
{
  ELFF::ClearErrPrefix();
}

} // extern "C"
