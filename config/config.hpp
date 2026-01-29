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


// Support out-of-source builds: if ELFF_CONFIG_FILE is defined, include it.
//
// Otherwise, use the local file: _config.hpp.

#ifndef ELFF_CONFIG_HPP
#define ELFF_CONFIG_HPP

#ifdef ELFF_CONFIG_FILE
#include ELFF_CONFIG_FILE
#else
#include "_config.hpp"
#endif

namespace ELFF
{

#if (defined(ELFF_USE_CUDA) && defined(__CUDACC__)) || \
    (defined(ELFF_USE_HIP) && defined(__HIPCC__))
#define ELFF_HOST_DEVICE __host__ __device__
#else
#define ELFF_HOST_DEVICE
#endif

// BEAM precision configuration

#if defined ELFF_USE_SINGLE && defined ELFF_USE_DOUBLE
#error "DOUBLE and SINGLE precision cannot both be specified"
#endif

#ifdef ELFF_USE_SINGLE
typedef float real_t;
#elif defined ELFF_USE_DOUBLE
typedef double real_t;
#else
#error "Either DOUBLE or SINGLE precision must be specified"
#endif


} // namespace beam

// Return value for main function in examples that should be skipped by testing
// in some case. This return value prevents failures in testing.
#define ELFF_SKIP_RETURN_VALUE 242

// Request a global object to be instantiated for each thread in its TLS.
#define ELFF_THREAD_LOCAL thread_local

// ELFF_DEPRECATED macro to mark obsolete functions and methods
// see https://stackoverflow.com/questions/295120/c-mark-as-deprecated
#if defined(__GNUC__) || defined(__clang__)
#define ELFF_DEPRECATED __attribute__((deprecated))
#elif defined(_MSC_VER)
#define ELFF_DEPRECATED __declspec(deprecated)
#else
#pragma message("WARNING: You need to implement ELFF_DEPRECATED for this compiler")
#define ELFF_DEPRECATED
#endif

// Common configuration macros


// Macro ELFF_EXPORT: this macro is used when declaring exported global
// variables and static class variables in public header files, e.g.:
//    extern ELFF_EXPORT Geometry Geometries;
//    static ELFF_EXPORT Device device_singleton;
// In cases where a class contains multiple static variables, instead of marking
// all such variables with ELFF_EXPORT, one can mark the class with ELFF_EXPORT,
// e.g.:
//    class ELFF_EXPORT MemoryManager ...
// Note: ELFF's GitHub CI includes a shared MSVC build that will fail if a
// variable that needs ELFF_EXPORT does not have it. However, builds with
// optional external libraries are not tested and may require separate checks to
// determine the necessity of ELFF_EXPORT.
#if defined(_MSC_VER) && defined(ELFF_SHARED_BUILD)
#ifdef beam_EXPORTS
#define ELFF_EXPORT __declspec(dllexport)
#else
#define ELFF_EXPORT __declspec(dllimport)
#endif
#else
#define ELFF_EXPORT
#endif
// On Cygwin the option -std=c++11 prevents the definition of M_PI. Defining
// the following macro allows us to get M_PI and some needed functions, e.g.
// posix_memalign(), strdup(), strerror_r().
#ifdef __CYGWIN__
#define _XOPEN_SOURCE 600
#endif

#endif // ELFF_CONFIG_HPP
