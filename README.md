# Basilisk Beam Library

[![Unit Tests](https://github.com/tortotubus/beam/actions/workflows/unit-tests.yml/badge.svg)](https://github.com/tortotubus/beam/actions/workflows/unit-tests.yml)


## Dependencies
For RPM-based distros (e.g. Redhat, Fedora):
```bash
sudo dnf install -y cmake gcc gcc-c++ git gcc-gfortran make gawk gnuplot doxygen gtest-devel openmpi-devel eigen3-devel hdf5-devel hdf5-openmpi-devel
```

For apt-based distros (e.g. Ubuntu):
```bash
sudo apt-get install -y  cmake gcc g++ git make gawk gnuplot doxygen libgtest-dev libopenmpi-dev libeigen3-dev libhdf5-dev  
```

## Building
From the repository root:
```bash
cmake -S . --preset Default 
cmake --build --preset Default --parallel  
```
 

## Basilisk model coupling

`library/ibm` contains the distributed mesh-manager implementation and its
supporting headers directly. Existing solver includes keep their timestep, kernel,
periodic-geometry, and output behavior. Dense registration remains the default;
`ibmeshmanager_set_model_sparse()` opts into sparse node storage under MPI.

Force-coupled models receive `nforce` from each node's center owner and publish
`npos`, `nvel`, and `nweight`. Support replicas do not contribute duplicate force
values. Sparse adapters implement `input_node` and `sync_node`, keyed by stable
`node_lid`, so the model owner need not hold all (or any) local IB nodes. Set the
mesh depth before sparse registration. Registration and synchronization must be
called in the same mesh order on every rank in `MPI_COMM_WORLD`.

ELFF force adapters own their transfer buffers and borrow their physical model
from IBRuntime. Before transferring an adapter to the manager, use
`elff_force_coupled_handle(adapter)` to configure its beam or rigid body; the
adapter's `ctx` is no longer the physical model pointer. Destroy adapters before
destroying the runtime. Model-state migration between owners is a separate step;
this change supplies force-coupled field transfer and sparse adapter support.

With `ELFF_BUILD_TESTS` and Basilisk enabled, build `force_coupled` and
`force_coupled_mpi`, then run:

```bash
ctest --test-dir build/debug -R BasiliskForceCoupled --output-on-failure
```

These compare dense/sparse beam trajectories with direct model integration,
including a nonzero model owner, overlapping MPI replicas, and periodic boundaries.

MPI case loading uses `ib_runtime_create_from_case_file_mpi(path, comm, ...)`.
It duplicates the communicator before parsing and constructing models, and
reports local case/model failures to every participating rank. MPI must already
be initialized, and all members of `comm` must call the factory. The original
`ib_runtime_create_from_case_file()` remains the serial entry point. Basilisk's
`elff_init()` selects the appropriate factory automatically; manual model
construction through `elff_runtime_get()` also establishes the communicator
first and stops with a diagnostic if setup fails.

XML-created models default to physical construction on runtime rank zero.
`<Models distribution="rank0">` explicitly selects that default;
`<Models distribution="round-robin">` assigns model `i` in XML order to rank
`i % communicator_size`. This balances model counts without considering fluid
partitions or model cost. Serial construction assigns every model to rank zero.
Ownership is fixed at construction; there is no per-model `ownerRank` setting.
Every rank receives the ordered descriptor directory;
remote views have a null physical `handle` and a valid `descriptor`.
`elff_register_model()` uses that descriptor to register sparse Basilisk storage.
Refinement levels come from the Basilisk case configuration, not model XML.
Free the mesh manager before destroying the runtime.

Mixed-owner runtime checkpoints collect state in model order and write one
file on communicator rank zero. Restore loads each current owner's models, so
the checkpoint is independent of the construction-time ownership policy.
All runtime ranks must participate in checkpointing. The standalone runtime
PolyData writer still requires a single model owner; Basilisk's mesh output
continues to use its distributed coupling data.

This path does not migrate physical models. The `.ib` checkpoint format now
stores logical mesh/node order (version 2) and restores local sparse support from
saved geometry. Legacy version-1 dense files in registration order can also be
loaded. ELFF restore refreshes model adapter caches before restoring IB fields.
The flapping-fibre case can recover an older incomplete checkpoint missing
`dump.ib` by restoring `dump.elff` and reconstructing nodal coupling force from
the restored fluid. It reports an approximate restart because the original
nodal force was not saved. Corrupt or incompatible sidecars remain errors.
The `.dt` sidecar records the actual iteration and appends the CFL timestep
history before output-time alignment, so restarting preserves timestep growth.
Older `.dt` files remain readable, using their saved timestep as a fallback for
that missing history. Checkpoint messages report the saved and next scheduled
iteration/time; the last terminal progress line may precede the saved state.

The flapping-fibre example parses `<Basilisk>` without DTD validation. It still
requires the configuration values it uses and checks their types and numeric
constraints. Additional elements or attributes and reordered children are
accepted. IBRuntime continues to validate the `<ELFF>` subtree against its DTD.

### Settling-fibre XML cases

`single_settling_fibre_{serial,mpi}` and `many_settling_fibres_{serial,mpi}`
accept one XML case path, like `flapping_fibre`. Examples are in
`basilisk/examples/cases/single_settling_fibre.xml` and
`basilisk/examples/cases/many_settling_fibres.xml`. From the repository root:

```sh
mpirun -np 2 build/release/basilisk/examples/single_settling_fibre_mpi basilisk/examples/cases/single_settling_fibre.xml
mpirun -np 2 build/release/basilisk/examples/many_settling_fibres_mpi basilisk/examples/cases/many_settling_fibres.xml
```

The single case requires exactly one model and retains its quiescent tank
boundaries. The many case uses a triply periodic box and one explicit
`EulerBeamHuang` entry per fibre; the example contains four fibres. Material
values in these examples are computed from the original Banaei parameters.
`Models distribution="round-robin"` enables static model balancing; omission
keeps all physical models on rank zero. Model geometry and material belong in
`ELFF`; domain, fluid, time, grid, uniform force density and output settings
belong in `Basilisk`. Both use the same Basilisk parser as flapping-fibre.

The XML path uses runtime-owned sparse models and writes MPI-reduced positions
and velocities to `fibre-centres.csv`. Output intervals are simulation-time
intervals. The existing CLI parameter paths, including random cloud generation,
remain available when no XML path is supplied. Settling coupling retains its
existing defaults (relaxation 0.8, two Richardson iterations, implicit bending).

Short fixtures in `tests/data/*settling*_short.xml` run to time 0.0003 with nine
nodes per fibre and reduced grid cost. Their CTest cases exercise serial and
two-rank MPI execution, with round-robin ownership in the many-fibre fixture.

### Case parameters and in-memory rendering

The flapping-fibre example in `basilisk/examples/cases/flapping_fibre.xml`
(and its template copy) uses this path. Edit `Re`, `massRatio`,
`bendingStiffness`, and `streamwiseForceDensity` in `Parameters`; the defaults
remain 200, 1.5, 0.001, and 0.75. Its unit fibre length, inlet speed, and fluid
density give `dynamicViscosity = 1/Re`. Numerical controls such as `beamNodes`,
`endTime`, and grid levels are also named parameters. The initial direction
retains its existing numeric components.

Cases can define inputs and derived coefficients as direct children of `Case`:

```xml
<Parameters>
  <Parameter name="Re" value="200"/>
  <Parameter name="beamNodes" expression="3*3"/>
</Parameters>
<Nondimensionalization>
  <Coefficient name="viscosity" expression="1/Re"/>
</Nondimensionalization>
```

Use these values in downstream attributes, for example
`<Geometry length="1" nodes="${parameters.beamNodes}"/>` or
`<Fluid density="1" dynamicViscosity="${coefficients.viscosity}">`.
Bare names such as `${beamNodes}` also refer to parameters. An expression must
occupy the whole attribute; string interpolation is not supported. Expressions
support arithmetic, right-associative `^`, `sqrt`, and `pi`. Parameters may
reference later parameters, but coefficients reference parameters only. Values
must use consistent units; no automatic unit conversion is performed.

The flapping and settling XML entry points render once per rank in memory and
pass the same rendered document to the Basilisk configuration reader and the
ELFF model factory. Rank zero writes `case.rendered.xml` alongside `case.xml`
for inspection. Neither reader loads that diagnostic file; failure to write it
only produces a warning. The original file is unchanged. The diagnostic retains
the parameter definitions and replaces downstream expression attributes with
round-trip-precision numbers. Invalid expressions stop construction with the
source line and attribute in the diagnostic.

Library callers can use `elff_xml_case_render` with a loaded document, or
`elff_xml_case_load_rendered_file` with a diagnostic-output flag. Both return
an independently owned document, released with `elff_xml_document_destroy`.
Pass it to `ib_runtime_create_from_case_document` (or its `_mpi` variant) and
other document readers before releasing it. The runtime copies the information
it needs during construction. The existing runtime file factories also render
in memory before validating the ELFF subtree. Low-level XML loading remains a
raw parse, and legacy numeric cases need no parameter sections.
