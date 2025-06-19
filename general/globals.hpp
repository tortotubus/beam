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

#pragma once

#include "../config/config.hpp"
#include <iostream>

namespace beam
{

/// Simple extension of std::ostream.
/** This class adds the ability to enable and disable the stream. The associated
    std::streambuf and tied std::ostream can be replaced with that of any
    std::ostream. */
class OutStream : public std::ostream
{
protected:
   // Pointer that stores the associated streambuf when output is disabled.
   std::streambuf *m_rdbuf;
   // Pointer that stores the tied ostream when output is disabled.
   std::ostream *m_tie;

   void Init();

public:
   /** @brief Construct an OutStream from the given stream @a os, by using its
       `rdbuf()`. */
   OutStream(std::ostream &os) : std::ostream(NULL) { SetStream(os); }

   /** @brief Replace the `rdbuf()` and `tie()` of the OutStream with that of
       @a os, enabling output. */
   void SetStream(std::ostream &os)
   {
      rdbuf(m_rdbuf = os.rdbuf()); tie(m_tie = os.tie()); Init();
   }

   /// Enable output.
   void Enable() { if (!IsEnabled()) { rdbuf(m_rdbuf); tie(m_tie); } }
   /// Disable output.
   void Disable()
   {
      if (IsEnabled()) { m_rdbuf = rdbuf(NULL); m_tie = tie(NULL); }
   }
   /// Check if output is enabled.
   bool IsEnabled() const { return (rdbuf() != NULL); }
};


/** @brief Global stream used by the library for standard output. Initially it
    uses the same std::streambuf as std::cout, however that can be changed.
    @sa OutStream. */
extern BEAM_EXPORT OutStream out;

/** @brief Global stream used by the library for standard error output.
    Initially it uses the same std::streambuf as std::cerr, however that can be
    changed.
    @sa OutStream. */
extern BEAM_EXPORT OutStream err;


/** @brief Construct a string of the form "<prefix><myid><suffix>" where the
    integer @a myid is padded with leading zeros to be at least @a width digits
    long. */
/** This is a convenience function, e.g. to redirect beam::out to individual
    files for each rank, one can use:
    \code
       std::ofstream out_file(MakeParFilename("app_out.", myid).c_str());
       beam::out.SetStream(out_file);
    \endcode
*/
std::string MakeParFilename(const std::string &prefix, const int myid,
                            const std::string suffix = "", const int width = 6);


/// @brief Wrapper for std::getenv.
///
/// @note Directly calling getenv causes a warning with MSVC. Use this wrapper
/// to suppress the warning.
const char* GetEnv(const char* name);

} // namespace beam
