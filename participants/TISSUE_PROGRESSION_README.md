# Progressive tissue-modelling RONS cases (Sections 4.1–4.4)

## Shared geometry, exposure, and agents

- Canonical cylinder: `inputs/cylinder/cylinder.msh`
- Radius `0.5 mm`, height `1.0 mm`, axis `+z`
- Treated `top` patch: H₂O₂ = NO₂⁻ = `0.1 mol/m3` (`100 µM`)
- Initial H₂O₂ = NO₂⁻ = `0 mol/m3`
- Cases 4.1–4.3: impermeable `bottom` and `walls` (`zeroGradient`)
- Case 4.4: absorbing tissue-domain `bottom` and `walls` (`fixedValue 0`)
- 5000 identical 20 µm packed agents in every case
- Spheroid center `(0, 0, 0.5 mm)`; agent-center z-range `0.317–0.683 mm`

The OpenFOAM scalar dimensions are molar concentration:
`[0 -3 0 0 1 0 0]`. BioDynaMo stores the same SI values.

## Transport coefficients

Both species deliberately use the same `D` and `k` within each run.

| Case | Environment | D for H₂O₂ and NO₂⁻ | k for H₂O₂ and NO₂⁻ | Duration |
|------|-------------|----------------------|----------------------|----------|
| **4.1** | in vitro | `1.75e-9 m2/s` active; `1.50e-11` and `1.50e-9` commented alternatives | `0 s-1` | `180 s` |
| **4.2** | in vitro | `1.75e-9 m2/s` | `0.005 s-1` active; `0.001` and `0.025` commented alternatives | `180 s` |
| **4.3** | in vitro | `1.75e-9 m2/s` | `0.005 s-1` | `180 s` |
| **4.4** | in vivo, matched to 4.3 | `1.75e-9 m2/s` | `0.005 s-1` | `900 s` |

The case-4.2 decay candidates correspond to half-lives of `693.1 s`
(`k=0.001`), `138.6 s` (`k=0.005`), and `27.7 s` (`k=0.025`). They are
sensitivity values: a single equal decay
constant for H₂O₂ and NO₂⁻ is a controlled-comparison assumption, not a claim
that the two species have identical reaction chemistry.

Cases 4.3 and 4.4 use the same provisional first-order uptake coefficient for
both species, `alpha_ext=0.02 s-1`. Each agent integrates separate
intracellular H₂O₂ and NO₂⁻ loads with `beta_R=0.02 s-1`, then integrates
damage with `alpha_D=2.5 m3/(mol s)` and `beta_D=0.02 s-1`. These
coefficients are uncalibrated placeholders pending the D3.2 values.
`alpha_D` is raised so that, at the intracellular loads reached during a
`180 s` in-vitro exposure, the more exposed agents cross `D=1`. Crossing
that threshold sets phenotype 2 but does not delete the agent. Each accepted
step then removes it with probability `1 - exp(-0.02 dt)`.

## Coupling progression

OpenFOAM solves:

```
∂Yi/∂t − ∇·(Di ∇Yi) + ki Yi = −YiSink
```

- **4.1:** diffusion only, one-way serial-explicit.
- **4.2:** extracellular decay, one-way serial-explicit.
- **4.3:** intracellular RONS loading and damage, one-way; uptake is not
  returned to OpenFOAM. Agents with damage at or above `1.0` are removed
  probabilistically after each accepted step.
- **4.4:** the same transport, intracellular ODEs, uptake rates, and removal
  hazard, with two-way H₂O₂ and NO₂⁻ feedback through `H2O2Sink` and
  `NO2Sink`. The draw is applied only after a preCICE window is accepted, so
  implicit-coupling rollback does not consume it. The in-vivo duration is
  still `900 s`, so this hazard removes a much larger fraction than in case 4.3.

OpenFOAM and the agent time step are `0.1 s`. Cases 4.1 and 4.2 use one
preCICE window of `180 s`, equal to the run duration, because the agents only
sample the extracellular fields. Cases 4.3 and 4.4 keep coupling windows of
`0.1 s`. OpenFOAM and agent visualization fields are written every `5 s`.

## How to run

Install OpenFOAM, preCICE, and BioDynaMo, then build and run a case, using the
commands in the repository [`README.md`](../README.md). The Make targets
are `run-tissue_41_diffusion`, `run-tissue_42_decay`,
`run-tissue_43_cell_response`, and `run-tissue_44_twoway`.

## Result figures (suggested)

- **4.1** Spatial H₂O₂ / NO₂⁻ fields for the three diffusivity values
- **4.2** Depth attenuation for the three decay values
- **4.3** Agent stress/damage/fate vs depth and time (extracellular field unchanged by cells)
- **4.4** Extracellular fields with local H₂O₂/NO₂⁻ cellular sinks vs 4.3
