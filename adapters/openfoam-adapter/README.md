# OpenFOAM–preCICE adapter: tissue subset

This is a reduced derivative of the preCICE OpenFOAM adapter for the four
tissue-modelling report examples. The adapter core is retained for preCICE
participant setup, volume meshes, boundary meshes, data exchange, and implicit
checkpointing. Only the `FA` module is built.

The reduced FA module:

- writes `H2O2` and `NO2` from finite-volume cell centres;
- reads `H2O2Sink` and `NO2Sink` into scalar finite-volume fields;
- writes the static boundary face and triangle data needed by the cells
  participant.

The adapter built here exchanges H2O2, NO2, H2O2Sink, NO2Sink, and the
boundary geometry used by the cells participant.

Build after loading the OpenFOAM v2406 environment:

```bash
./Allwmake
```

The upstream license is retained in `LICENSE`, and upstream citation metadata
is retained in `CITATION.cff`.
