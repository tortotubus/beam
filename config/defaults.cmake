if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
      "Build type: Debug, Release, RelWithDebInfo, or MinSizeRel." FORCE)
endif()

# # BEAM options. Set to mimic the default "defaults.mk" file.
option(BUILD_SHARED_LIBS "Enable shared library build of BEAM" OFF)
# option(BEAM_USE_MPI "Enable MPI parallel build" OFF)
# option(BEAM_USE_METIS "Enable METIS usage" ${BEAM_USE_MPI})
set(BEAM_PRECISION "double" CACHE STRING
    "Floating-point precision to use: single, or double")
option(BEAM_USE_EXCEPTIONS "Enable the use of exceptions" OFF)
# option(BEAM_USE_ZLIB "Enable zlib for compressed data streams." OFF)
# option(BEAM_USE_LIBUNWIND "Enable backtrace for errors." OFF)
# option(BEAM_USE_LAPACK "Enable LAPACK usage" OFF)
# option(BEAM_THREAD_SAFE "Enable thread safety" OFF)
# option(BEAM_USE_OPENMP "Enable the OpenMP backend" OFF)
# option(BEAM_USE_LEGACY_OPENMP "Enable legacy OpenMP usage" OFF)
# option(BEAM_USE_MEMALLOC "Enable the internal MEMALLOC option." ON)
# option(BEAM_USE_SUNDIALS "Enable SUNDIALS usage" OFF)
# option(BEAM_USE_SUITESPARSE "Enable SuiteSparse usage" OFF)
# option(BEAM_USE_SUPERLU "Enable SuperLU_DIST usage" OFF)
# option(BEAM_USE_SUPERLU5 "Use the old SuperLU_DIST 5.1 version" OFF)
# option(BEAM_USE_MUMPS "Enable MUMPS usage" OFF)
# option(BEAM_USE_STRUMPACK "Enable STRUMPACK usage" OFF)
# option(BEAM_USE_GINKGO "Enable Ginkgo usage" OFF)
# option(BEAM_USE_AMGX "Enable AmgX usage" OFF)
# option(BEAM_USE_MAGMA "Enable MAGMA usage" OFF)
# option(BEAM_USE_GNUTLS "Enable GNUTLS usage" OFF)
# option(BEAM_USE_GSLIB "Enable GSLIB usage" OFF)
# option(BEAM_USE_HDF5 "Enable HDF5 usage" OFF)
# option(BEAM_USE_NETCDF "Enable NETCDF usage" OFF)
# option(BEAM_USE_PETSC "Enable PETSc support." OFF)
# option(BEAM_USE_SLEPC "Enable SLEPc support." OFF)
# option(BEAM_USE_MPFR "Enable MPFR usage." OFF)
# option(BEAM_USE_SIDRE "Enable Axom/Sidre usage" OFF)
# option(BEAM_USE_FMS "Enable FMS usage" OFF)
# option(BEAM_USE_CONDUIT "Enable Conduit usage" OFF)
# option(BEAM_USE_PUMI "Enable PUMI" OFF)
# option(BEAM_USE_HIOP "Enable HiOp" OFF)
# option(BEAM_USE_CUDA "Enable CUDA" OFF)
# option(BEAM_USE_HIP "Enable HIP" OFF)
# option(BEAM_USE_OCCA "Enable OCCA" OFF)
# option(BEAM_USE_RAJA "Enable RAJA" OFF)
# option(BEAM_USE_CEED "Enable CEED" OFF)
# option(BEAM_USE_UMPIRE "Enable Umpire" OFF)
# option(BEAM_USE_SIMD "Enable use of SIMD intrinsics" OFF)
# option(BEAM_USE_ADIOS2 "Enable ADIOS2" OFF)
# option(BEAM_USE_CALIPER "Enable Caliper support" OFF)
# option(BEAM_USE_ALGOIM "Enable Algoim support" OFF)
# option(BEAM_USE_MKL_CPARDISO "Enable MKL CPardiso" OFF)
# option(BEAM_USE_MKL_PARDISO "Enable MKL Pardiso" OFF)
# option(BEAM_USE_ADFORWARD "Enable forward mode for AD" OFF)
# option(BEAM_USE_CODIPACK "Enable automatic differentiation (AD) using CoDiPack" OFF)
# option(BEAM_USE_BENCHMARK "Enable Google Benchmark" OFF)
# option(BEAM_USE_PARELAG "Enable ParELAG" OFF)
# option(BEAM_USE_TRIBOL "Enable Tribol" OFF)
# option(BEAM_USE_ENZYME "Enable Enzyme" OFF)
# # Optional overrides for autodetected MPIEXEC and MPIEXEC_NUMPROC_FLAG
# # set(BEAM_MPIEXEC "mpirun" CACHE STRING "Command for running MPI tests")
# # set(BEAM_MPIEXEC_NP "-np" CACHE STRING
# #     "Flag for setting the number of MPI tasks")
# set(BEAM_MPI_NP 4 CACHE STRING "Number of processes used for MPI tests")
# # Allow a user to disable testing, examples, and/or miniapps at CONFIGURE TIME
# # if they don't want/need them (e.g. if BEAM is "just a dependency" and all they
# # need is the library, building all that stuff adds unnecessary overhead). Note
# # that the examples or miniapps can always be built using the targets 'examples'
# # or 'miniapps', respectively.
# option(BEAM_ENABLE_TESTING "Enable the ctest framework for testing" ON)
# option(BEAM_ENABLE_EXAMPLES "Build all of the examples" OFF)
# #option(BEAM_ENABLE_MINIAPPS "Build all of the miniapps" OFF)
# option(BEAM_ENABLE_BENCHMARKS "Build all of the benchmarks" OFF)
# # Setting CXX/MPICXX on the command line or in user.cmake will overwrite the
# # autodetected C++ compiler.
# # set(CXX g++)
# # set(MPICXX mpicxx)
# # Set the target CUDA architecture
# set(CUDA_ARCH "sm_60" CACHE STRING "Target CUDA architecture.")
# # Set the target HIP architecture
# set(HIP_ARCH "gfx900" CACHE STRING "Target HIP architecture.")
# set(BEAM_DIR ${CMAKE_CURRENT_SOURCE_DIR})
# # The *_DIR paths below will be the first place searched for the corresponding
# # headers and library. If these fail, then standard cmake search is performed.
# # Note: if the variables are already in the cache, they are not overwritten.
# set(HYPRE_DIR "${BEAM_DIR}/../hypre/src/hypre" CACHE PATH
#     "Path to the hypre library.")
# # If hypre was compiled to depend on BLAS and LAPACK:
# # set(HYPRE_REQUIRED_PACKAGES "BLAS" "LAPACK" CACHE STRING
# #     "Packages that HYPRE depends on.")
# # CUDA and HIP dependencies for HYPRE are handled in FindHYPRE.cmake.
# set(METIS_DIR "${BEAM_DIR}/../metis-4.0" CACHE PATH "Path to the METIS library.")
# set(LIBUNWIND_DIR "" CACHE PATH "Path to Libunwind.")
# # For sundials_nvecmpiplusx and nvecparallel remember to build with MPI_ENABLE=ON
# # and modify cmake variables for hypre for sundials
# set(SUNDIALS_DIR "${BEAM_DIR}/../sundials-5.0.0/instdir" CACHE PATH
#     "Path to the SUNDIALS library.")
# # The following may be necessary, if SUNDIALS was built with KLU:
# # set(SUNDIALS_REQUIRED_PACKAGES "SuiteSparse/KLU/AMD/BTF/COLAMD/config"
# #     CACHE STRING "Additional packages required by SUNDIALS.")
# set(SuiteSparse_DIR "${BEAM_DIR}/../SuiteSparse" CACHE PATH
#     "Path to the SuiteSparse library.")
# set(SuiteSparse_REQUIRED_PACKAGES "BLAS" "METIS"
#     CACHE STRING "Additional packages required by SuiteSparse.")
# set(ParMETIS_DIR "${BEAM_DIR}/../parmetis-4.0.3" CACHE PATH
#     "Path to the ParMETIS library.")
# set(ParMETIS_REQUIRED_PACKAGES "METIS" CACHE STRING
#     "Additional packages required by ParMETIS.")
# set(SuperLUDist_DIR "${BEAM_DIR}/../SuperLU_DIST_8.1.2" CACHE PATH
#     "Path to the SuperLU_DIST library.")
# # SuperLU_DIST may also depend on "OpenMP", depending on how it was compiled.
# set(SuperLUDist_REQUIRED_PACKAGES "MPI" "ParMETIS" "METIS"
#     "LAPACK" "BLAS" CACHE STRING
#     "Additional packages required by SuperLU_DIST.")
# set(MUMPS_DIR "${BEAM_DIR}/../MUMPS_5.5.0" CACHE PATH
#     "Path to the MUMPS library.")
# # MUMPS may also depend on "OpenMP", depending on how it was compiled.
# set(MUMPS_REQUIRED_PACKAGES "MPI" "MPI_Fortran" "ParMETIS" "METIS"
#     "ScaLAPACK" "LAPACK" "BLAS" CACHE STRING
#     "Additional packages required by MUMPS.")
# # If the MPI package does not find all required Fortran libraries:
# # set(MUMPS_REQUIRED_LIBRARIES "gfortran" "mpi_mpifh" CACHE STRING
# #     "Additional libraries required by MUMPS.")
# set(STRUMPACK_DIR "${BEAM_DIR}/../STRUMPACK-build" CACHE PATH
#     "Path to the STRUMPACK library.")
# # STRUMPACK may also depend on "OpenMP", depending on how it was compiled.
# # Starting with v2.2.0 of STRUMPACK, ParMETIS and Scotch are optional.
# set(STRUMPACK_REQUIRED_PACKAGES "MPI" "MPI_Fortran" "ParMETIS" "METIS"
#     "Scotch/ptscotch/ptscotcherr/scotch/scotcherr"
#     "ScaLAPACK" "LAPACK" "BLAS" CACHE STRING
#     "Additional packages required by STRUMPACK.")
# # If the MPI package does not find all required Fortran libraries:
# # set(STRUMPACK_REQUIRED_LIBRARIES "gfortran" "mpi_mpifh" CACHE STRING
# #     "Additional libraries required by STRUMPACK.")
# # The Scotch library, required by STRUMPACK <= v2.1.0, optional in STRUMPACK >=
# # v2.2.0.
# set(Scotch_DIR "${BEAM_DIR}/../scotch_6.0.4" CACHE PATH
#     "Path to the Scotch and PT-Scotch libraries.")
# set(Scotch_REQUIRED_PACKAGES "Threads" CACHE STRING
#     "Additional packages required by Scotch.")
# # Tell the "Threads" package/module to prefer pthreads.
# set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
# set(Threads_LIB_VARS CMAKE_THREAD_LIBS_INIT)
# # The ScaLAPACK library, required by STRUMPACK
# set(ScaLAPACK_DIR "${BEAM_DIR}/../scalapack-2.0.2/lib/cmake/scalapack-2.0.2"
#     CACHE PATH "Path to the configuration file scalapack-config.cmake")
# set(ScaLAPACK_TARGET_NAMES scalapack)
# # set(ScaLAPACK_TARGET_FORCE)
# # set(ScaLAPACK_IMPORT_CONFIG DEBUG)
# set(Ginkgo_DIR "${BEAM_DIR}/../ginkgo" CACHE PATH "Path to the Ginkgo library.")
# set(AMGX_DIR "${BEAM_DIR}/../amgx" CACHE PATH "Path to AmgX")
# set(MAGMA_DIR "${BEAM_DIR}/../magma" CACHE PATH "Path to MAGMA")
# set(MAGMA_REQUIRED_PACKAGES "BLAS" "LAPACK" CACHE STRING
#     "Additional packages required by MAGMA.")
# set(GNUTLS_DIR "" CACHE PATH "Path to the GnuTLS library.")
# set(GSLIB_DIR "" CACHE PATH "Path to the GSLIB library.")
# set(HDF5_DIR "/usr" CACHE PATH "Path to the HDF5 library.")
# set(NETCDF_DIR "" CACHE PATH "Path to the NetCDF library.")
# set(NetCDF_REQUIRED_PACKAGES "HDF5/C/HL" CACHE STRING
#     "Additional packages required by NetCDF.")
# set(PETSC_DIR "${BEAM_DIR}/../petsc" CACHE PATH
#     "Path to the PETSc main directory.")
# set(PETSC_ARCH "arch-linux2-c-debug" CACHE STRING "PETSc build architecture.")
# set(SLEPC_DIR "${BEAM_DIR}/../slepc" CACHE PATH
#     "Path to the SLEPc main directory.")
# set(SLEPC_ARCH "arch-linux2-c-debug" CACHE STRING "SLEPC build architecture.")
# set(MPFR_DIR "" CACHE PATH "Path to the MPFR library.")
# set(FMS_DIR "${BEAM_DIR}/../fms" CACHE PATH
#     "Path to the FMS library.")
# # If FMS is built with Conduit:
# # set(FMS_REQUIRED_PACKAGES "Conduit/relay" CACHE STRING
# #     "Additional packages required by FMS.")
# set(CONDUIT_DIR "${BEAM_DIR}/../conduit" CACHE PATH
#     "Path to the Conduit library.")
# set(AXOM_DIR "${BEAM_DIR}/../axom" CACHE PATH "Path to the Axom library.")
# # May need to add "Boost" as requirement.
# if (BEAM_USE_SIDRE)
#     if (BEAM_USE_MPI)
#         set(Axom_REQUIRED_PACKAGES "Conduit/blueprint/blueprint_mpi/relay/relay_mpi" CACHE STRING
#             "Additional packages required by Axom.")
#     elseif()
#         set(Axom_REQUIRED_PACKAGES "Conduit/blueprint/relay" CACHE STRING
#             "Additional packages required by Axom.")
#     endif()
# endif()
# set(PUMI_DIR "${BEAM_DIR}/../pumi-2.1.0" CACHE STRING
#     "Directory where PUMI is installed")
# set(HIOP_DIR "${BEAM_DIR}/../hiop/install" CACHE STRING
#     "Directory where HiOp is installed")
# set(HIOP_REQUIRED_PACKAGES "BLAS" "LAPACK" CACHE STRING
#     "Packages that HiOp depends on.")
# set(MKL_CPARDISO_DIR "" CACHE STRING "MKL installation path.")
# set(MKL_MPI_WRAPPER_LIB "mkl_blacs_mpich_lp64" CACHE STRING "MKL MPI wrapper library")
# set(MKL_LIBRARY_DIR "" CACHE STRING "Custom library subdirectory")
# set(MKL_PARDISO_DIR "" CACHE STRING "MKL installation path.")
# set(OCCA_DIR "${BEAM_DIR}/../occa" CACHE PATH "Path to OCCA")
# set(RAJA_DIR "${BEAM_DIR}/../raja" CACHE PATH "Path to RAJA")
# set(CEED_DIR "${BEAM_DIR}/../libCEED" CACHE PATH "Path to libCEED")
# set(UMPIRE_DIR "${BEAM_DIR}/../umpire" CACHE PATH "Path to Umpire")
# set(CALIPER_DIR "${BEAM_DIR}/../caliper" CACHE PATH "Path to Caliper")
# set(BLITZ_DIR "${BEAM_DIR}/../blitz" CACHE PATH "Path to Blitz")
# set(ALGOIM_DIR "${BEAM_DIR}/../algoim" CACHE PATH "Path to Algoim")
# set(Algoim_REQUIRED_PACKAGES "Blitz" CACHE STRING
#     "Packages that ALGOIM depends on.")
# set(BENCHMARK_DIR "${BEAM_DIR}/../google-benchmark" CACHE PATH
#     "Path to Google Benchmark")
# # Provide paths, since ParELAG is dependent on BEAM and BEAM needs to be
# # compiled (or at least cmake needs to succeed) before compiling ParELAG.
# set(PARELAG_DIR "${BEAM_DIR}/../parelag" CACHE PATH "Path to ParELAG")
# set(PARELAG_INCLUDE_DIRS "${PARELAG_DIR}/src;${PARELAG_DIR}/build/src" CACHE
#     STRING "Path to ParELAG headers.")
# set(PARELAG_LIBRARIES "${PARELAG_DIR}/build/src/libParELAG.a" CACHE STRING
#     "The ParELAG library.")
# set(TRIBOL_DIR "${BEAM_DIR}/../tribol" CACHE PATH "Path to Tribol")
# set(Tribol_REQUIRED_PACKAGES "Axom/core/mint/slam/slic" CACHE STRING
#     "Additional packages required by Tribol")
# set(ENZYME_DIR "${BEAM_DIR}/../enzyme" CACHE PATH "Path to Enzyme")
# set(BLAS_INCLUDE_DIRS "" CACHE STRING "Path to BLAS headers.")
# set(BLAS_LIBRARIES "" CACHE STRING "The BLAS library.")
# set(LAPACK_INCLUDE_DIRS "" CACHE STRING "Path to LAPACK headers.")
# set(LAPACK_LIBRARIES "" CACHE STRING "The LAPACK library.")
# set(CODIPACK_INCLUDE_DIRS "${BEAM_DIR}/../CoDiPack/include" CACHE STRING
#     "Path to CoDiPack headers.")
# set(CODIPACK_LIBRARIES "")

# # Some useful variables:
# # set(CMAKE_SKIP_PREPROCESSED_SOURCE_RULES ON) # Skip *.i rules
# # set(CMAKE_SKIP_ASSEMBLY_SOURCE_RULES  ON)    # Skip *.s rules
# # set(CMAKE_VERBOSE_MAKEFILE ON CACHE BOOL "Verbose makefiles.")
