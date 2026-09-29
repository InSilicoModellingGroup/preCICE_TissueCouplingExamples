# Public report subset

This repository is intentionally narrower than the development repository.
Its purpose is to accompany the four tissue-modelling examples described in
report Sections 4.1–4.4 without publishing the separate benchmark feature set.

## Allowed coupling data

The OpenFOAM FA module may expose only:

- volume scalars written by OpenFOAM: `H2O2`, `NO2`;
- volume scalar sinks read by OpenFOAM: `H2O2Sink`, `NO2Sink`;
- stationary boundary geometry used by the cells participant:
  `WallGeometry`, `WallNormals`, `WallAreas`, `BoundaryTriangleVertex0`,
  `BoundaryTriangleVertex1`, `BoundaryTriangleVertex2`, `BoundaryKind`,
  `BoundaryNormals`, and `BoundaryVelocity`.

`BoundaryVelocity` is attached to the exact boundary-triangle representation;
for these static tissue cases it is zero. It is not a carrier-flow coupling
field.

## Explicitly excluded

- particle-to-fluid momentum or semi-implicit carrier feedback in the public
  coupling interface;
- particle mass feedback in the public coupling interface;
- OpenFOAM FA exchange of fluid velocity, pressure, acceleration, vorticity,
  or turbulence fields;
- CHT, FSI, and fluid-fluid adapter modules;
- cavity, venturi, aneurysm, oxygen-cube, and other benchmark participants;
- generated OpenFOAM meshes, processor directories, time results, BioDynaMo
  outputs, logs, profiling data, local dependencies, and report source files.

The scalar sink feedback in case 4.4 is retained because it is the defining
coupling step of the report example.

The cells participant and BDM adapter retain shared lifecycle, mapping, and
cell-model infrastructure needed to compile and run these four examples. Their
active preCICE data contract is limited to the fields listed above.
