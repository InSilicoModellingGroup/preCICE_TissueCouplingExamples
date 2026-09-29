# BioDynaMo–preCICE adapter: tissue subset

This directory contains the header-only adapter portion required by the four
tissue examples. It maps the OpenFOAM volume fields `H2O2` and `NO2` to
BioDynaMo agents, imports static boundary geometry, and maps cellular uptake
back to the `H2O2Sink` and `NO2Sink` fields for case 4.4.

Particle momentum and particle mass feedback are outside the report scope and
have been removed from this repository.

The cells application includes `src/bdm_adapter.h` directly. Supporting
headers provide parameter structures, mesh mapping, cell-interface traits, and
boundary-triangle reconstruction.
