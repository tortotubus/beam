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

#ifdef _WIN32
// Turn off CRT deprecation warnings for getenv
#define _CRT_SECURE_NO_WARNINGS
#endif

#include "../config/config.hpp"
#include "globals.hpp"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cstdlib>  // getenv

namespace beam
{

OutStream out(std::cout);
OutStream err(std::cerr);

namespace internal
{
bool beam_out_initialized = false;
bool beam_err_initialized = false;
}

void OutStream::Init()
{
   if (this == &beam::out)
   {
      internal::beam_out_initialized = true;
   }
   else if (this == &beam::err)
   {
      internal::beam_err_initialized = true;
   }
}

std::string MakeParFilename(const std::string &prefix, const int myid,
                            const std::string suffix, const int width)
{
   std::stringstream fname;
   fname << prefix << std::setw(width) << std::setfill('0') << myid << suffix;
   return fname.str();
}


#ifdef BEAM_USE_MPI

MPI_Comm BEAM_COMM_WORLD = MPI_COMM_WORLD;

MPI_Comm GetGlobalMPI_Comm()
{
   return BEAM_COMM_WORLD;
}

void SetGlobalMPI_Comm(MPI_Comm comm)
{
   BEAM_COMM_WORLD = comm;
}

#endif

const char *GetEnv(const char* name)
{
   return getenv(name);
}

}
