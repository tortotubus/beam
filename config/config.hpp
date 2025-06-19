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


// Support out-of-source builds: if BEAM_CONFIG_FILE is defined, include it.
//
// Otherwise, use the local file: _config.hpp.

#ifndef BEAM_CONFIG_HPP
#define BEAM_CONFIG_HPP

#ifdef BEAM_CONFIG_FILE
#include BEAM_CONFIG_FILE
#else
#include "_config.hpp"
#endif

namespace beam
{

#if (defined(BEAM_USE_CUDA) && defined(__CUDACC__)) || \
    (defined(BEAM_USE_HIP) && defined(__HIPCC__))
#define BEAM_HOST_DEVICE __host__ __device__
#else
#define BEAM_HOST_DEVICE
#endif

// BEAM precision configuration

#if defined BEAM_USE_SINGLE && defined BEAM_USE_DOUBLE
#error "DOUBLE and SINGLE precision cannot both be specified"
#endif

#ifdef BEAM_USE_SINGLE
typedef float real_t;
#elif defined BEAM_USE_DOUBLE
typedef double real_t;
#else
#error "Either DOUBLE or SINGLE precision must be specified"
#endif


} // namespace beam

// Return value for main function in examples that should be skipped by testing
// in some case. This return value prevents failures in testing.
#define BEAM_SKIP_RETURN_VALUE 242

// Request a global object to be instantiated for each thread in its TLS.
#define BEAM_THREAD_LOCAL thread_local

// BEAM_DEPRECATED macro to mark obsolete functions and methods
// see https://stackoverflow.com/questions/295120/c-mark-as-deprecated
#if defined(__GNUC__) || defined(__clang__)
#define BEAM_DEPRECATED __attribute__((deprecated))
#elif defined(_MSC_VER)
#define BEAM_DEPRECATED __declspec(deprecated)
#else
#pragma message("WARNING: You need to implement BEAM_DEPRECATED for this compiler")
#define BEAM_DEPRECATED
#endif

// Common configuration macros


// Macro BEAM_EXPORT: this macro is used when declaring exported global
// variables and static class variables in public header files, e.g.:
//    extern BEAM_EXPORT Geometry Geometries;
//    static BEAM_EXPORT Device device_singleton;
// In cases where a class contains multiple static variables, instead of marking
// all such variables with BEAM_EXPORT, one can mark the class with BEAM_EXPORT,
// e.g.:
//    class BEAM_EXPORT MemoryManager ...
// Note: BEAM's GitHub CI includes a shared MSVC build that will fail if a
// variable that needs BEAM_EXPORT does not have it. However, builds with
// optional external libraries are not tested and may require separate checks to
// determine the necessity of BEAM_EXPORT.
#if defined(_MSC_VER) && defined(BEAM_SHARED_BUILD)
#ifdef beam_EXPORTS
#define BEAM_EXPORT __declspec(dllexport)
#else
#define BEAM_EXPORT __declspec(dllimport)
#endif
#else
#define BEAM_EXPORT
#endif
// On Cygwin the option -std=c++11 prevents the definition of M_PI. Defining
// the following macro allows us to get M_PI and some needed functions, e.g.
// posix_memalign(), strdup(), strerror_r().
#ifdef __CYGWIN__
#define _XOPEN_SOURCE 600
#endif

#endif // BEAM_CONFIG_HPP
