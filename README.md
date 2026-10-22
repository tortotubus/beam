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
from ModelRuntime. Before transferring an adapter to the manager, use
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

MPI case loading uses `model_runtime_create_from_case_file_mpi(path, comm, ...)`.
It duplicates the communicator before parsing and constructing models, and
reports local case/model failures to every participating rank. MPI must already
be initialized, and all members of `comm` must call the factory. The original
`model_runtime_create_from_case_file()` remains the serial entry point. Basilisk's
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
Coupled restarts require both `dump.elff` and `dump.ib` alongside the fluid
`dump`. All ranks check required files before restoring any state; missing,
corrupt, or incompatible required files fail restart. Missing IB state is not
reconstructed approximately.
The `.dt` sidecar records the actual iteration and appends the CFL timestep
history before output-time alignment, so restarting preserves timestep growth.
Older `.dt` files remain readable, using their saved timestep as a fallback for
that missing history. Checkpoint messages report the saved and next scheduled
iteration/time; the last terminal progress line may precede the saved state.

The flapping-fibre executables now require the v2 DTD, `caseType="flapping_fibre"`
and `Grid type="quadtree"`. Their shared reader in `io/xml/xml.h` loads
boundaries, adaptation, coupling, restart, and output settings. The thin
case parser reads the optional `inletSpeed` parameter (default one) for
fresh-start velocity; inlet boundary values remain explicit XML settings.
The case retains tip statistics and computes `omega` before fluid output.

Output field names are bound after solver initialization; unknown names fail
with an XML source location. Fluid and solid schedules default to an interval
of one and the run end, with the first output at one interval. Omitted output
or restart blocks are disabled. `Checkpoint` retains signal/wall-clock defaults;
an optional `SimulationTime` adds the existing checkpoint service's timestep
trigger. Additional checkpoint-field selection is currently rejected explicitly.
Kernel names are `Cosine` (default) and `Roma`. The IB level can exceed the fluid
adaptation ceiling, as in the original flapping configuration.

### Settling-fibre XML cases

The `settling_fibre_{serial,mpi}` drivers share the v2 reader and events in
`basilisk/library/io/xml/xml.h`. Both cases require triply periodic boundaries
and accept any fibre model supported by the ELFF XML factory; fluid velocity
starts at zero on a fresh run. Model configuration belongs in the ELFF XML block.

```sh
cmake --build build/debug --target settling_fibre_serial settling_fibre_mpi
build/debug/basilisk/examples/settling_fibre_serial basilisk/examples/settling_fibre_1.xml
mpirun -np 2 build/debug/basilisk/examples/settling_fibre_mpi basilisk/examples/settling_fibre_1000.xml
```

The single example has one fixed centred horizontal fibre. The many example
requests 1000 randomly positioned horizontal fibres. Change `count`,
`anchorMin*`/`anchorMax*`, and `Models distribution` to control insertion and
MPI ownership. Clearance checks are nonperiodic: keep initial geometry away
from periodic faces if clearance across those faces matters.

Both XML files retain the Banaei expressions: `A = pi*length^2/(4*rp^2)`,
`W = r*rho0*A*g`, linear density `(1+r)*rho0*A`, bending rigidity
`gamma*r*W*length^3`, dynamic viscosity `rho0*sqrt(r*g*length^3)/Ga`, and
force density `(0,-W,0)`. Both currently use `gamma=0.01`. These formulas
remain in XML, including the original extra factor `r` in bending rigidity.

`Grid/AdaptWavelet` supplies the velocity tolerance; `ImmersedBoundary/ForceCoupled`
supplies `Iterations` (2 in the examples) and `Relaxation` (0.8). Named
`Parameters` supply `statsInterval` (a positive integer, default 1), which
controls iteration-based reporting
of MPI-reduced centres and velocities to `fibre-centres.csv`.

The v2 `Output` blocks select fluid and solid fields and their simulation-time
intervals. `l2` and `omega_z` are refreshed before fluid output. `Restart` enables
restoration; a checkpoint requires the complete `dump`, `dump.elff`, and
`dump.ib` set. `Checkpoint` enables the shared checkpoint handler.

### Case parameters and in-memory rendering

The flapping-fibre example in `basilisk/examples/cases/flapping_fibre.xml`
(and its template copy) uses this path. Edit `Re`, `massRatio`,
`bendingStiffness`, and `streamwiseForceDensity` in `Parameters`; the defaults
remain 200, 1.5, 0.001, and 0.75. Its unit fibre length, inlet speed, and fluid
density give `dynamicViscosity = 1/Re`. Numerical controls such as `beamNodes`,
`endTime`, and grid levels are also named parameters. The initial direction
retains its existing numeric components.

Cases can define inputs and derived parameters as direct children of `Case`:

```xml
<Parameters>
  <Parameter name="Re" value="200"/>
  <Parameter name="beamNodes" expression="3*3"/>
  <Parameter name="viscosity" expression="1/Re"/>
</Parameters>
```

Use these values in downstream attributes, for example
`<Geometry length="1" nodes="${parameters.beamNodes}"/>` or
`<DynamicViscosity>${parameters.viscosity}</DynamicViscosity>`.
Bare names such as `${beamNodes}` also refer to parameters. An expression must
occupy the whole attribute; string interpolation is not supported. Expressions
support arithmetic, right-associative `^`, `sqrt`, and `pi`. Parameters may
reference later parameters; cyclic references are rejected. Values
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
Pass it to `model_runtime_create_from_case_document` (or its `_mpi` variant) and
other document readers before releasing it. The runtime copies the information
it needs during construction. The existing runtime file factories also render
in memory before validating the ELFF subtree. Low-level XML loading remains a
raw parse, and legacy numeric cases need no parameter sections.

### Ownership of Basilisk case XML code

`basilisk/library/io/input-case.h` is the reusable, header-only reader for
common `<Basilisk>` settings. Its `BasiliskCaseConfiguration` and
`basilisk_case_*` functions cover domain, time, fluid, grid, output, and optional
immersed-boundary level/force settings. They do not interpret `caseType`,
construct models, use simulation globals, or call MPI. A solver without any
fibre or ELFF model can use the same reader.

Each template has exactly one case-specific parser header:
`flapping_fibre_parser.h`, `single_settling_fibre_parser.h`, or
`many_settling_fibres_parser.h`. Types and functions use that case's name.
Each case configuration contains a `BasiliskCaseConfiguration basilisk` member;
this composition makes the shared settings explicit. These headers initialize
the ELFF runtime and enforce the case's model-count requirements: flapping and
many-settling require at least one model; single-settling requires exactly one.
They are included after the template's Basilisk/ELFF runtime headers.

`basilisk_case_read` borrows an already rendered document. `basilisk_case_load`
loads and renders once, reads the common settings, and returns an owned document
on success (null on failure). The case initializer agrees on parsing success
across MPI ranks, passes that same document to ELFF, and then destroys it.
Only the designated writer requests the optional `.rendered.xml` diagnostic.

Fibre-specific registration, force-refresh events, and fibre CSV statistics
belong to each `.c` template. The two settling templates intentionally keep
local copies rather than depending on another case or adding fibre-specific
helpers to the library. Generic expression evaluation and XML document ownership
remain in ELFF's XML library. Case parser headers never include each other.

### Centralized XML schemas

`io/DTD/ELFF/ELFF-v1.dtd` contains the complete authoritative schema, with no
external includes. Both models reuse `Geometry`, `Material`, `InitialLine`,
`Anchor`, and `Direction`; each element is declared exactly once.

`InitialLine` contains an `Anchor` and a `Direction`. Huang requires
`anchorEnd="left|right"`. The implicit rod accepts omission or `anchorEnd="left"`:
its anchor always identifies node zero, and direction points toward the last
node. These model-specific restrictions are checked by the semantic parsers;
the shared DTD declares `anchorEnd` optional. Rod XML now uses `InitialLine`
in place of `RodInitialLine`.

`InitialLine` also supports explicit `RandomAnchor` and `RandomDirection` choices
for insertion groups. See [XML placement](docs/insertion.md#xml-placement) for
`InsertionDomain`, group counts, retry budgets, and serial/MPI examples.

The XML element is `<EulerBeamHuang>`, matching the numerical class.
Its construction helper is `EulerBeamHuangBuilder.hpp`, the IB
wrapper is `IBEulerBeamHuang`, and C/Basilisk constructors use
`ib_euler_beam_huang_*` / `elff_euler_beam_huang_*`.
Beam solver class names omit `Inextensible` consistently (for example,
`EulerBeamHuang`, `EulerBeamGGL`, and `EulerBeamPenalty`). This naming change
does not alter the models’ physical constraints.

Edit `io/DTD/ELFF/ELFF-v1.dtd` for XML contracts. To add another supported
model, add its declarations and update the Models alternatives in that file. There are no per-model CMake registrations. Runtime readers and constructors
still need to support those model tags. Huang's XML reader lives in
`io/case/models/EulerBeamHuangXML.hpp/.cpp`. Its `parse` method returns an owned
`EulerBeamHuangDefinition`; `parse_placement` reads fixed or permitted random
placement separately. Neither method constructs models or assigns MPI ranks.
The returned values remain valid after the XML document is destroyed.
`models/beam/EulerBeamHuangBuilder` builds the initial mesh and boundary values
from that definition and a resolved `InitialLineDefinition`. These data types
and the builder can also be used without XML. Shared owned placement types live
in `io/case/ModelPlacement.hpp` under `ELFF::IO::Case`; XML placement parsing lives in `io/case/InitialLineXML.hpp` under
`ELFF::IO::Case`.
The implicit rod reader likewise lives in `io/case/models/BendingRodImplicitXML.hpp`,
with owned settings in `models/articulated/BendingRodImplicitDefinition.hpp`.
Its builder depends on those settings, without including the XML parser.

`io/CXX/xmlCase.hpp` provides a small C++ view over the C XML API: `XMLNode`
borrows a node and provides typed attribute/child reads; `XMLError` preserves
the C status and source diagnostic. The document must outlive its node views.
The wrapper does not load, render, or own documents. Model-specific validation
remains in each model parser. The construction layer in
`c/models/ModelRuntimeCase.cpp` dispatches to these parsers, adds model indices to
errors, and constructs models in XML order. `ModelRuntime` owns the resulting
models; its `IBCoupling` component registers their coupling views. See
[Runtime ownership and IB coupling](docs/manual/runtime.md) for the C and C++ APIs.

At configure time, `config/scripts/bundle_dtd.py` checks the schema and writes
`build/<configuration>/generated/ELFF-v1.dtd`. Python 3.7 or newer is required,
with no third-party Python packages. CMake tracks the source and script so schema
edits regenerate the artifact on the next build. The runtime embeds this
artifact, and installation copies it to `share/elff/dtd/ELFF/ELFF-v1.dtd`.
The source, embedded, and installed schemas need no external entity resolver.

DTD element names are global: shared contracts must be declared once, while
incompatible model-specific contracts require distinct names. The centralized DTD now covers the complete
`Case`: optional `Parameters`, followed by
`ELFF`, `Basilisk`, or both (in either order). It also supports standalone
section validation. Basilisk declares the common fluid/grid/output settings,
optional immersed-boundary settings, and optional settling controls.

Call `elff_xml_case_validate(document, &error)` to validate the complete case
against the embedded schema, before or after expression rendering. Existing
section readers retain their scope: the runtime factory validates ELFF, and the
low-level Basilisk reader performs semantic checks and permits extensions.
The flapping and settling templates validate the complete rendered case before
collective model construction. Whole-case validation rejects undeclared extensions; add their declarations to
the centralized DTD when introducing a new case configuration.

For editor completion, a case in `basilisk/examples` can declare:

```xml
<!DOCTYPE Case SYSTEM "../../io/DTD/ELFF/ELFF-v1.dtd">
```

Numeric attributes use `CDATA` so source expressions remain valid. Positivity,
numeric ranges, parameter value/expression exclusivity, and case-type-specific
requirements remain semantic checks. Schema versions remain 1.

### ELFF version-2 XML

ELFF's case factories accept `<Case formatVersion="2" dimension="3">`
(or `dimension="2"`). The dimension is required; the `ELFF` section has no
version or dimension attributes. Version 1 remains available. Both DTDs are
embedded and installed.

Version 2 uses `InitialCondition/Line` and numeric material children:
`Material/BendingRigidity` (force × length²) and
`Material/LinearDensity` (mass/length). Numeric attributes and leaf text
accept whole-value `${expression}` substitutions. Parameter
definitions remain available after rendering.

Each model has two boundary elements with distinct `end="left|right"`
attributes: `EulerBeamBoundary` accepts `free|pinned|clamped`, while
`BendingRodBoundary` supports a free right end and a pinned or free left end.
`BendingRodSegments` selects the implicit rod. Huang still requires
`Line anchorEnd`; rods permit only left anchoring or omission.

[The v2 rod insertion example](examples/cases/insert_rods_v2.xml) works with
the XML insertion driver. Model ordering, MPI ownership, insertion behavior,
and checkpoint payloads are unchanged.

### Implicit rod physical damping

In v2 rod XML, append optional `BendingDamping` to `Material`:

```xml
<Material>
  <BendingRigidity>0.01</BendingRigidity>
  <LinearDensity>1</LinearDensity>
  <BendingDamping>0.001</BendingDamping>
</Material>
```

The coefficient defaults to zero and must be finite and nonnegative. In SI its
units are N m² s; each joint receives rotational damping `D/h`, where `h` is
segment length. This uses the existing projected relative-angular-velocity
joint damping, with equal and opposite torques and nonpositive damping power.
It damps internal bending motion without damping uniform translation. It is
separate from generalized-alpha numerical damping (`spectralRadius`) and is
not a remedy for every partitioned fluid-coupling instability. Huang does not
support this element. C++ callers set `damping_rigidity`; explicit per-joint
`bending_damping` values override it. Checkpoints retain the resolved joint
coefficients in the existing payload format.

### Standalone stochastic isotropic flow

`basilisk/templates/isotropic_fibre.c` is now a **fluid-only** experiment using
Basilisk's vanilla `navier-stokes/centered.h`. Despite the historical filename,
it has no fibre, ELFF, MDF, XML, or custom checkpoint dependency. Both isotropic
drivers use `grid/octree.h`. The small plain-C forcing helper lives beside the
template in `isotropic_forcing.h`.

```sh
cmake --build build/debug --target isotropic_fibre_serial isotropic_fibre_mpi
# Run in an output directory: energy.csv is overwritten on each fresh run.
/path/to/isotropic_fibre_serial 4 10 0.05 1 42 5
```

Arguments are `minlevel end_time sigma correlation_time seed maxlevel output_interval`, with
trailing arguments optional. Defaults are `4 10 0.05 1 1 4 0.05`. Equal grid levels
disable adaptation. The cube has side `2*pi`, viscosity `0.01`, maximum timestep
`0.001`, and adaptation tolerance `0.01`; these constants can be edited in the
template. The fluid starts at rest and the OU process starts stationary.
HTG output uses the existing HDF5 writer, so build these drivers with the CMake
targets above (which supply the HDF5/MPI headers and libraries). This output
helper does not introduce ELFF model or immersed-boundary solver dependencies.
`isotropic_flow_output/htg.vtkhdf.series` contains pressure `p`, velocity `u`,
and the lambda-2 vortex diagnostic `l2`, starting at time zero and then every
`output_interval` simulation-time units. Open the series file in ParaView.

Nine conjugate pairs cover integer wavevectors with `0 < |n|^2 <= 2`. Rank zero
updates all Gaussian OU coefficients with the actual fluid timestep and
broadcasts them. Transverse projection and direct face-centre evaluation give
continuous periodic, divergence-free acceleration without FFTs. The sampled
field is added to `a`; its separate face field is retained for diagnostics.

With variance `sigma^2` for **each real and imaginary** coefficient component
and a factor of two for conjugate reconstruction, the expected spatial mean of
squared force is `72*sigma^2`, giving total RMS `6*sqrt(2)*sigma`. A formula of
`36*sigma^2` instead assumes variance `sigma^2/2` for each real component.
No mode normalization or constant-power adjustment is applied.

`energy.csv` records end-of-step time, actual dt, total specific kinetic energy,
midpoint external forcing work, strain-based viscous dissipation, sampled force
RMS, and discrete force-divergence RMS. All averages use cell volumes and MPI
reductions. The work pairs velocity across the final correction with the
metric-weighted face-to-cell forcing; it is a diagnostic, not an exact discrete
budget closure including pressure, advection, viscosity and adaptation errors.
No IBM power exists in this fluid-only case. The seed reproduces a realization
for a given timestep sequence; changing adaptive timesteps changes that sequence.
No restart support or statistically stationary turbulence calibration is included.

### Per-case Basilisk defaults and XML overrides

Flapping and settling now initialize a complete case configuration before reading
v2 XML. An ELFF-only document is accepted; `<Basilisk/>` has the same effect as
omitting the block. Model selection and structural parameters still belong in
`<ELFF>` and remain required.

| Default | Flapping | Settling |
|---|---|---|
| Domain | length 16, origin (-2,-8) | length 32, origin (-16,-16,-16) |
| Fluid density / viscosity | 1 / 0.005 | 1 / sqrt(0.1)/40 |
| Maximum dt / end time | 0.0003 / 50 | 0.001 / 10 |
| Initial/minimum / maximum / IB level | 6 / 10 / 11 | 4 / 12 / 12 |
| Boundaries | unit x inlet and pressure outlet | periodic on all three axes |
| Uniform force density | (0.75,0,0) | (0,-0.1*pi/(4*30²),0) |
| Velocity adaptation tolerance | 0.005 | 0.001*sqrt(0.1) |
| Force iterations / relaxation | 2 / 1 | 2 / 0.8 |
| Fluid/solid output interval | 0.05 | 0.1 |
| Output directory | flapping_fibre_output | settling_fibre_output |

Both default to the cosine kernel, enabled restart and checkpoint handling, and
their previous fluid/solid field selections. A missing checkpoint starts a fresh
run; an incomplete checkpoint set still fails. Defaults are simulation-unit
values, not automatically recomputed from named parameters or model material.
Use XML expressions when parameters should determine fluid settings.

Only explicitly supplied values override the baseline, for example:

```xml
<Basilisk>
  <Time end="2"/>
  <Fluid>
    <FluidMaterial><DynamicViscosity>0.02</DynamicViscosity></FluidMaterial>
  </Fluid>
  <Grid><Depth>4</Depth><AdaptWavelet enable="false"/></Grid>
  <Restart enable="false"/>
  <ImmersedBoundary level="5"/>
  <Output directory="trial">
    <FluidOutput><SimulationTime step="0.1"/></FluidOutput>
    <SolidOutput enable="false"/>
    <Checkpoint enable="false"/>
  </Output>
</Basilisk>
```

A supplied boundary-condition list replaces all default boundary entries and
periodic axes; a Fluid block containing only material settings preserves them.
Explicit adaptation or output field entries replace the corresponding default
field list. A schedule-only output block preserves its fields. Omitted blocks
inherit their enabled state; `enable="false"` explicitly disables adaptation,
output, restart, or checkpoint handling. Disabling adaptation sets the grid's
minimum and maximum to its initial depth. With adaptation enabled, an overridden
Depth must lie within the merged refinement range. Output schedule ends follow
the merged run duration unless explicitly supplied.

Optional `caseType` and grid `type` remain executable assertions when present.
The final merged configuration is validated before model construction; invalid
explicit values fail instead of silently falling back. FluidMaterial accepts
Density/DynamicViscosity child elements or the existing density/dynamicViscosity
attributes, but rejects specifying the same property both ways.

### XML case launcher

On Linux with XML support enabled, `elff` selects the flapping or settling driver
from the required `Basilisk.caseType` attribute. It locates the driver beside
its own executable; grid type is checked by the selected driver.

```sh
elff settling_fibre_1.xml
mpiexec -n 8 elff settling_fibre_1.xml
srun -n 8 elff settling_fibre_1.xml
```

Serial/MPI selection is automatic from Open MPI/PMI/PMIx rank metadata or a
Slurm job-step rank, including single-rank launches. A Slurm allocation alone
does not select MPI. The dispatcher never starts a launcher or initializes MPI.
Inherited launcher metadata can also trigger MPI selection.

Arguments after the XML path are passed unchanged to the driver. There are no
case or execution-mode overrides; files omitting `Basilisk` can still be passed
directly to a case executable. Unsupported case types and missing executables
produce errors. The fluid-only isotropic example uses positional arguments
rather than this XML interface.
