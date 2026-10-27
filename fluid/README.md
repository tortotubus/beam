# Fluid solver

The WIP incompressible AMR solver builds independently of ELFF. The external
Grains adapter, its runtime options, and its stepping/load-transfer hooks have
been removed. Standalone DLM sphere constraints and hydrodynamic loads remain.
Public headers live under `elff/fluid/` and use `ELFF::Fluid`, with `DLM` and
`Boundary` subnamespaces. Internal helpers use `SolverDetail` and `DLM::Detail`. Solver parameters use the `fluid` ParmParse prefix
(for example, `fluid.fixed_dt`).

The `pacific` multiplier layout is a numerical discretization choice and does
not require the Grains library.

Build using the checked-out AMReX and AMReX-Hydro submodules:

```sh
cmake -S fluid -B build/fluid -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build/fluid --target elff_fluid_solver -j 8
build/fluid/elff_fluid_solver fluid/inputs/taylor-green
```

Defaults are 3D, double precision, CPU, serial, without embedded boundaries.
AMReX cache options can enable MPI/OpenMP or select precision. Configure a fresh
build directory when changing dependency configurations. To use installed
packages instead, set `ELFF_FLUID_BUNDLED_AMREX=OFF` and supply their prefixes
through `CMAKE_PREFIX_PATH`.

Verified with the checked-out AMReX `008b67375045` and Hydro `154572aa6cf2`:
the library and driver compile/link, and two-step Taylor–Green and zero-flow
prescribed-sphere DLM runs complete in the default configuration. These smoke
runs do not establish numerical accuracy or validate MPI/GPU builds.

The checked-out AMReX requires CMake 3.25 and publishes `cxx_std_20` on its
targets. The solver requests C++17; linking AMReX raises the effective standard
to C++20. This does not require raising the standard for unrelated ELFF targets.

## ELFF integration work

- Add an optional fluid component to the root build and link its implementation
  into ELFF, or export a separate `ELFF::fluid` library. Currently the root source
  collection does not include fluid. Keep the executable separate from the
  library and decide whether fluid is a separately installed component.
- Extend root header copying and installation to include `elff/fluid/*.H`: the
  root installer currently selects `.h` and `.hpp`, so it would omit these `.H`
  files. The public API uses the `ELFF::Fluid` namespace.
- Align `amrex::Real` with ELFF's configurable `real_t`, preferably with a
  compile-time type check at the boundary. Match SINGLE/DOUBLE in both dependency
  builds to `ELFF_PRECISION`; do not silently convert internal solution fields.
- Define AMReX initialization/finalization and MPI communicator ownership for an
  embedded solver. The standalone executable owns AMReX's lifetime; an ELFF API
  must accommodate an application that already initialized MPI or AMReX.
- Connect configuration to ELFF's case/runtime model if the solver should run
  through the ELFF launcher. Currently input is AMReX ParmParse, output is AMReX
  plotfiles, and body geometry is standalone spheres. ELFF model/IBM coupling,
  XML contracts, and any HDF5/checkpoint integration require explicit adapters.
- Add numerical regressions to the existing GoogleTest executables: projection
  divergence, viscous decay, AMR synchronization/regridding, and standalone DLM
  motion/loads. Compile and exercise supported precision, dimension, and parallel
  configurations before promising those configurations.

## AMReX external-project work

The root external-project wrappers are not used by this standalone build and
still need integration:

- Correct Hydro's submodule source to
  `third_party/amrex_hydro/amrex_hydro_submodule` and give it distinct download
  and binary directories. It currently uses AMReX's source and build paths.
- Order Hydro after AMReX's install (`DEPENDS amrex_ep` when both are built), and
  pass AMReX's install prefix/package location to Hydro's configure step.
- Forward toolchains/compilers, build type, install libdir, PIC/shared settings,
  precision, dimension, MPI, OpenMP, GPU backend, and relevant feature options.
  Hydro needs projections; AMReX needs AMR and linear solvers. EB is optional for
  this solver, although Hydro enables it by default. Fortran and particles are
  not needed. Keep matching Hydro/AMReX options consistent.
- Provide imported targets and build dependencies to consuming ELFF targets,
  including include directories and actual library byproducts. An
  `ExternalProject_Add` target alone is not a linkable CMake library target.
  Adding to `CMAKE_PREFIX_PATH` inside a child directory also does not propagate
  package discovery to the parent. `find_package` cannot discover an install
  that will only be produced later during the same build.
- Either use a true two-stage superbuild with `find_package` in the inner build,
  or configure both dependencies with `add_subdirectory` and consume their real
  targets, as this standalone build does. The latter avoids manually recreating
  transitive compiler/link requirements. Pin compatible revisions instead of
  downloading moving `development` branches.
- Support already installed dependencies when third-party builds are disabled,
  enforce Hydro's dependency on AMReX, and add dependency discovery to the
  installed ELFF package config if exported targets expose them. Check a clean
  install with a separate consumer; root ELFF currently neither links these
  libraries nor declares them in `ELFFConfig.cmake.in`.
