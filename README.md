# preCICE tissue-coupling examples

This repository is the report-facing subset of the FVM–ABM coupling work. It
contains the four progressive tissue-modelling examples, the BioDynaMo `cells`
participant, the scalar diffusion solver, and only the adapter paths needed by
these examples.

The larger development repository is intentionally not mirrored here. In
particular, this subset does not contain the other flow benchmarks or particle
momentum feedback.

## Examples

| Report case | Directory | Coupling |
|---|---|---|
| 4.1 diffusion | `tissue_41_diffusion` | H2O2/NO2 transport, one-way |
| 4.2 decay | `tissue_42_decay` | transport plus extracellular decay, one-way |
| 4.3 cell response | `tissue_43_cell_response` | intracellular RONS response, one-way |
| 4.4 two-way | `tissue_44_twoway` | intracellular response plus H2O2/NO2 sink feedback |

All cases use the same cylindrical domain and initial 5000-cell spheroid. See
[`modular_Coupling/participants/TISSUE_PROGRESSION_README.md`](modular_Coupling/participants/TISSUE_PROGRESSION_README.md)
for the equations, coefficients, and progression between cases.

## Included code

- `applications/myDiffusionFoam`: transient multi-species diffusion/decay
  solver with optional scalar sink fields.
- `adapters/openfoam-adapter`: OpenFOAM adapter core and the reduced FA module.
  Its report data contract is H2O2, NO2, H2O2Sink, NO2Sink, and boundary
  geometry only.
- `adapters/bdm-adapter`: header-only preCICE bridge used by the cells
  participant for scalar sampling, boundary geometry, and cellular sinks.
- `modular_Coupling/cells`: BioDynaMo participant and the four tissue parameter
  sets/initial agent layouts.
- `modular_Coupling/participants/tissue_*`: clean OpenFOAM inputs and portable
  preCICE configuration templates.

The exact inclusion/exclusion policy is recorded in [`SCOPE.md`](SCOPE.md).

## Requirements

- OpenFOAM v2406
- preCICE 3.x with development headers and `pkg-config` metadata
- BioDynaMo v1.05.143 (or a compatible installation)
- CMake 3.19+, MPI, Gmsh, and Python 3

Dependencies and generated simulation results are not vendored.

## Build

```bash
source /usr/lib/openfoam/openfoam2406/etc/bashrc

cd adapters/openfoam-adapter
./Allwmake

cd ../../applications/myDiffusionFoam
wmake

cd ../../modular_Coupling/cells
make build BDM_PATH=/path/to/biodynamo/bin/thisbdm.sh
```

## Run a case

Prepare and start OpenFOAM in the first terminal:

```bash
cd modular_Coupling/participants/tissue_41_diffusion
./config-coupling.sh
./run.sh
```

After OpenFOAM starts waiting for its coupling partner, start cells in a second
terminal:

```bash
cd modular_Coupling/cells
make run-tissue_41_diffusion BDM_PATH=/path/to/biodynamo/bin/thisbdm.sh
```

Replace the case name with `tissue_42_decay`, `tissue_43_cell_response`, or
`tissue_44_twoway` as needed. `config-coupling.sh` generates an ignored runtime
preCICE XML file from the tracked `.xml.in` template, inserting an absolute
socket exchange directory without committing a machine-specific path.

## Publication check

Run this before committing or publishing:

```bash
./scripts/check-public-scope.sh
```

The check rejects excluded adapter modules, momentum-feedback identifiers,
absolute home-directory paths, generated time directories, and unexpected
participant cases.

## Licensing and attribution

The reduced OpenFOAM adapter derives from the preCICE OpenFOAM adapter and
retains its license and citation metadata in `adapters/openfoam-adapter/`.
Before making the repository public, add the intended top-level license for the
project-owned solver, BioDynaMo adapter, and example material.
