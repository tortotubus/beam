# Huang fibre in two-dimensional forced periodic flow

`isotropic_fibre.c` uses a quadtree, the Huang beam model, and the existing
`centered-mdf.h` immersed-boundary solver. The serial and MPI wrappers are in
`basilisk/examples`.

The forcing follows the local Basilisk example
`third_party/basilisk/basilisk_submodule/src/examples/isotropic.c`:
`a = forcing_rate * (u - mean(u))`, with volume-weighted mean velocity.
The default periodic square has side `2*pi`, viscosity `0.01`, and forcing
rate `0.1`. A deterministic sum of divergence-free Fourier modes replaces
the three-dimensional ABC initial condition. `U0` is the initial continuum
RMS speed. The straight fibre starts at rest, centred in the square, with
both ends free, zero gravity, and implicit Huang bending. Its initial
orientation `b_theta` is in radians; `huang_gamma` and `huang_rho` are passed
as bending stiffness and mass per unit length. These are configurable fibre
properties; the flapping-flow setup and its gravity are not used.

The default `beam.ds_dx = 2` matches the flapping case:
`ds/dx = (1/64)/(16/2^11) = 2`. With `beam.b_nodes = 0`, the node count is
computed from `b_length`, `L_fluid`, and `ibmlevel`, rounding the number of
segments to the nearest integer. Defaults give 21 nodes and `ds/dx ≈ 2.037`.
An explicit positive `b_nodes` overrides automatic spacing. The actual ratio
is printed at startup. `beam.force_relaxation` controls the existing IBM
Richardson update (default 0.1).

This is a 2D analogue, not a reproduction of the 3D isotropic benchmark.
The oriented fibre breaks statistical isotropy. Linear forcing in 2D can
feed box-scale structures; no large-scale drag or energy controller is
included, and statistically stationary turbulence is not guaranteed.
Refinement, timestep, domain size, and fibre resolution need convergence
checks for a physical study.

## Build and run

From the repository root with the usual build environment:

```bash
cmake -S . --preset Default
cmake --build --preset Default --target isotropic_fibre_serial isotropic_fibre_mpi
build/release/basilisk/examples/isotropic_fibre_serial --template > isotropic.params
build/release/basilisk/examples/isotropic_fibre_serial isotropic.params
```

For MPI, generate a separate parameter file (the default output directory
includes `mpi`) and run:

```bash
build/release/basilisk/examples/isotropic_fibre_mpi --template > isotropic-mpi.params
mpirun -np 2 build/release/basilisk/examples/isotropic_fibre_mpi isotropic-mpi.params
```

A short, coarse check is provided:

```bash
build/release/basilisk/examples/isotropic_fibre_serial basilisk/templates/isotropic_fibre_smoke.params
```

Use a fresh `base_path` for independent runs. The existing checkpoint/restore
handlers are enabled; CSV statistics append on restart. `experiment.t_end`
sets the stopping time. `basilisk.fluid.maxlevel` controls flow refinement;
`ibmlevel` independently controls refinement near the fibre. Levels are
restricted to 2–12 in this starter case.

## Output

- `flow.csv`: iteration, time, mean fluctuating kinetic energy, viscous
  dissipation, enstrophy, nominal continuum forcing power `2*forcing_rate*KE`,
  and both mean velocity components. Dissipation excludes beam energy and
  immersed-boundary work, so these columns alone are not a total energy budget.
- HDF/VTK flow output: vorticity, pressure, velocity, and immersed-boundary force.
- HDF/VTK fibre output: geometry, interpolated fluid velocity, force, and fibre
  velocity.

`output_interval` controls field output and `stats_interval` controls the CSV
sampling in timesteps. Defaults stop at time 10; the smoke case stops at 0.004.
