# preCICE tissue-coupling examples

Coupled simulations of extracellular reactive-oxygen and nitrogen species in a
cylindrical tissue domain. OpenFOAM (finite-volume) transports H2O2 and NO2.
BioDynaMo (agent-based) represents a 5000-cell spheroid. [preCICE](https://precice.org/)
exchanges the fields between the two programs.

The four cases below use the same geometry and the same initial spheroid. Each
case adds one modelling step to the previous one.

These instructions assume a fresh Ubuntu 22.04 machine. The commands use the
install locations produced by the steps in this file. Run the repository
commands from the root of this clone.

## 1. System packages

```bash
sudo apt-get update
sudo apt-get install -y \
  build-essential \
  ca-certificates \
  cmake \
  git \
  gmsh \
  pkg-config \
  python3 \
  wget \
  xvfb
```

`cmake` must be 3.19 or newer. Ubuntu 22.04 ships a new enough version.
`gmsh` converts the shared cylinder mesh the first time a case is configured.
`xvfb` lets the BioDynaMo participant write its ParaView state file on a
machine without a display.

## 2. preCICE 3.3.0

The OpenFOAM adapter and the cells participant both link against preCICE 3.x
and locate it with `pkg-config`.

```bash
wget https://github.com/precice/precice/releases/download/v3.3.0/libprecice3_3.3.0_jammy.deb
sudo apt-get install -y ./libprecice3_3.3.0_jammy.deb
pkg-config --modversion libprecice
```

`pkg-config` should print `3.3.0`.

## 3. OpenFOAM v2406

The case scripts load OpenFOAM from the official package location
`/usr/lib/openfoam/openfoam2406/etc/bashrc`.

```bash
wget -q -O - https://dl.openfoam.com/add-debian-repo.sh | sudo bash
sudo apt-get update
sudo apt-get install -y openfoam2406-dev
```

## 4. BioDynaMo

The cells participant is built against BioDynaMo commit
`a9d3c90e97164660d0ce567a357eb1cfe38035aa` (v1.05.143). The commands below
install it in `$HOME/biodynamo` and produce
`$HOME/biodynamo/build/bin/thisbdm.sh`.

`prerequisites.sh` uses `sudo` and installs BioDynaMo's own third-party
packages. It asks for confirmation unless `SILENT_INSTALL=1` is set.

```bash
git clone https://github.com/BioDynaMo/biodynamo.git "$HOME/biodynamo"
cd "$HOME/biodynamo"
git checkout a9d3c90e97164660d0ce567a357eb1cfe38035aa
SILENT_INSTALL=1 ./prerequisites.sh required
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j"$(nproc)"
```

`make` downloads and builds BioDynaMo's third-party libraries, including ROOT
and ParaView. The first build needs network access and takes a long time.

Check that the environment file exists:

```bash
test -f "$HOME/biodynamo/build/bin/thisbdm.sh" && echo "BioDynaMo environment file is ready"
```

Load it in every terminal that builds or runs the cells participant:

```bash
source "$HOME/biodynamo/build/bin/thisbdm.sh"
```

Sourcing prints a confirmation and sets `BDMSYS`.

## 5. Build the codes in this repository

The BioDynaMo adapter in `adapters/bdm-adapter` is header-only. It is compiled
as part of the cells participant. The OpenFOAM adapter and the diffusion
solver are built separately.

Open a terminal, go to the root of this repository, and run:

```bash
source /usr/lib/openfoam/openfoam2406/etc/bashrc

cd adapters/openfoam-adapter
./Allwmake

cd ../../applications/myDiffusionFoam
wmake
```

`Allwmake` writes `libpreciceAdapterFunctionObject.so` under
`adapters/openfoam-adapter/Make/linux64GccDPInt32Opt/`. The case scripts add
that directory to `LD_LIBRARY_PATH`. `wmake` installs `myDiffusionFoam` into
the OpenFOAM user application directory, which is on `PATH` after OpenFOAM is
loaded. The case scripts load OpenFOAM themselves, so that binary is found
when `./run.sh` starts.

In a second terminal, from the repository root:

```bash
source "$HOME/biodynamo/build/bin/thisbdm.sh"
cd participants/cells
make build
```

`make build` is enough to compile the participant. Repeat it after changing
the cells source. Each `make run-...` command also rebuilds before it starts
the simulation.

## 6. Run a coupled case

Each case uses two terminals. Start OpenFOAM first. It is the preCICE acceptor
and waits for the cells participant.

The shipped decompositions use 16 OpenFOAM ranks
(`numberOfSubdomains` in each case's `system/decomposeParDict`). On a machine
with fewer cores, set that value to at least 2, then run `./config-coupling.sh`
again before `./run.sh`.

### Terminal 1 — OpenFOAM

From the repository root, for the diffusion case:

```bash
cd participants/tissue_41_diffusion
./config-coupling.sh
./run.sh
```

`config-coupling.sh` imports the cylinder mesh on the first call, decomposes
the case, and writes the runtime preCICE configuration from the tracked
`.xml.in` template. `./run.sh` does not decompose the case. Run
`./config-coupling.sh` again after changing the rank count.

Wait until the OpenFOAM log prints:

```text
Setting up primary communication to coupling partner/s
```

### Terminal 2 — BioDynaMo

From the repository root:

```bash
source "$HOME/biodynamo/build/bin/thisbdm.sh"
cd participants/cells
make run-tissue_41_diffusion
```

Both participants then advance together. OpenFOAM writes
`participants/tissue_41_diffusion/logs/tissue_41_diffusion.log`
and, after the run, a reconstructed case that can be opened as
`tissue_41_diffusion.foam`. The cells participant writes
`participants/cells/output/tissue_41_diffusion/`.

Use the same two-terminal sequence for the other cases. Replace the directory
name and the Make target together:

| Case | Terminal 1 directory | Terminal 2 command |
|---|---|---|
| 4.1 diffusion | `participants/tissue_41_diffusion` | `make run-tissue_41_diffusion` |
| 4.2 decay | `participants/tissue_42_decay` | `make run-tissue_42_decay` |
| 4.3 cell response | `participants/tissue_43_cell_response` | `make run-tissue_43_cell_response` |
| 4.4 two-way | `participants/tissue_44_twoway` | `make run-tissue_44_twoway` |

Run one case at a time. The participants share
`precice-run/` for the preCICE connection files.

To remove generated OpenFOAM results from a case, run `./clean.sh` in that
case directory while the solver is stopped. `make clean` in
`participants/cells` removes the cells build and output files.

## The four cases

All four cases use the cylinder in `inputs/cylinder/`
(radius 0.5 mm, height 1.0 mm) and the same 5000-cell spheroid. The treated
top face is held at H2O2 = NO2 = 0.1 mol/m3. The interior starts at zero.
OpenFOAM advances with a time step of 0.1 s. Cases 4.1 and 4.2 use one
preCICE window of 180 s, the same as the run duration, because the spheroid
only reads the extracellular fields. Cases 4.3 and 4.4 keep 0.1 s windows so
the intracellular response, and in case 4.4 the sink feedback, stay in step
with the solver. Cases 4.1–4.3 run for 180 s. Case 4.4 runs for 900 s.

Coefficients, boundary values, and the intracellular equations are listed in
[`participants/TISSUE_PROGRESSION_README.md`](participants/TISSUE_PROGRESSION_README.md).

### 4.1 Diffusion — `tissue_41_diffusion`

OpenFOAM diffuses H2O2 and NO2 from the treated top face into the cylinder.
There is no extracellular decay and no cellular uptake. The spheroid only
reads the two concentration fields, so the case uses a single one-way
exchange over the full 180 s.

### 4.2 Decay — `tissue_42_decay`

The same diffusion problem, with an added first-order extracellular decay of
both species. The spheroid still only reads the fields, so this case also uses
a single one-way exchange over the full 180 s.

### 4.3 Cell response — `tissue_43_cell_response`

Transport and extracellular decay are unchanged. Each agent now integrates its
own intracellular H2O2 and NO2 loads, accumulates damage, and can be removed
after the damage threshold is crossed. That uptake is kept inside BioDynaMo
and is not sent back to OpenFOAM, so the extracellular fields do not depend on
the cells. The exchange stays one-way.

### 4.4 Two-way exchange — `tissue_44_twoway`

The intracellular model of case 4.3 is kept, and the cellular uptake rates are
written back to OpenFOAM as the volumetric sinks `H2O2Sink` and `NO2Sink`.
The side and bottom boundaries absorb the species. The run covers the longer
in-vivo duration of 900 s.

## Repository layout

```text
adapters/openfoam-adapter/     OpenFOAM–preCICE adapter (FA module)
adapters/bdm-adapter/          header-only BioDynaMo–preCICE adapter
applications/myDiffusionFoam/  multi-species diffusion/decay solver
participants/cells/            BioDynaMo participant and per-case parameters
participants/tissue_4*/        the four OpenFOAM cases
inputs/cylinder/               shared cylinder mesh
cpu-layout.sh                  core pinning used by the run scripts
```

## Licenses

`adapters/openfoam-adapter/` is derived from the preCICE OpenFOAM adapter and
keeps that project's license in `adapters/openfoam-adapter/LICENSE`. BioDynaMo
is Apache-2.0. preCICE and OpenFOAM keep their own licenses.
