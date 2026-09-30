#!/bin/bash

SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO_ROOT=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO_ROOT/cpu-layout.sh"

# ── Environment ───────────────────────────────────────────────────────────────
set +u
source /usr/lib/openfoam/openfoam2406/etc/bashrc
set -euo pipefail
export LD_LIBRARY_PATH="$REPO_ROOT/adapters/openfoam-adapter/Make/$WM_OPTIONS:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

cd "$SCRIPT_DIR"

# ── Prepared decomposition ────────────────────────────────────────────────────
NPROCS=$(awk '/numberOfSubdomains/{ gsub(";",""); print $2; exit }' system/decomposeParDict)
if [[ -z "$NPROCS" || "$NPROCS" -le 1 ]]; then
    echo "ERROR: tissue_42_decay coupling requires numberOfSubdomains > 1 in system/decomposeParDict" >&2
    exit 1
fi
if ! compgen -G "processor*" >/dev/null; then
    echo "ERROR: No processor directories found. Run ./config-coupling.sh first; run.sh does not decompose." >&2
    exit 1
fi
ACTUAL_NPROCS=$(find processor* -maxdepth 0 -type d 2>/dev/null | wc -l)
if [ "$ACTUAL_NPROCS" -ne "$NPROCS" ]; then
    echo "ERROR: Existing decomposition has $ACTUAL_NPROCS processors, but system/decomposeParDict requests $NPROCS." >&2
    echo "       Run ./config-coupling.sh after changing the rank count." >&2
    exit 1
fi
TISSUE_42_DECAY_CPUSET="${TISSUE_42_DECAY_CPUSET:-$(build_physical_cpuset_from_slice 0 "$NPROCS")}"
echo "→ Using prepared $NPROCS-way decomposition from config-coupling.sh"

# ── Coupled run ───────────────────────────────────────────────────────────────
mkdir -p logs
echo "→ Starting coupled run  (log → logs/tissue_42_decay.log)..."
echo "→ Pinning OpenFOAM ranks to physical cores: $TISSUE_42_DECAY_CPUSET"
if command -v taskset >/dev/null 2>&1; then
    taskset -c "$TISSUE_42_DECAY_CPUSET" mpirun --bind-to core --map-by core -np "$NPROCS" myDiffusionFoam -parallel 2>&1 | tee logs/tissue_42_decay.log
else
    echo "→ taskset not available; relying on mpirun core binding only"
    mpirun --bind-to core --map-by core -np "$NPROCS" myDiffusionFoam -parallel 2>&1 | tee logs/tissue_42_decay.log
fi
echo "→ Coupled run finished."

# ── Post-process ──────────────────────────────────────────────────────────────
HAS_WRITTEN_TIME=false
if [[ -d processor0 ]]; then
    for d in processor0/*; do
        [[ -d "$d" ]] || continue
        time_name=${d##*/}
        if [[ "$time_name" =~ ^[0-9]+([.][0-9]+)?$ ]] &&
           awk -v t="$time_name" 'BEGIN { exit !(t > 0) }'; then
            HAS_WRITTEN_TIME=true
            break
        fi
    done
fi

if [[ "$HAS_WRITTEN_TIME" == true ]]; then
    echo "→ Reconstructing case from processor directories..."
    reconstructPar 2>&1 | tee logs/tissue_42_decay.reconstruct.log
else
    echo "→ No nonzero processor time was written; skipping reconstruction."
fi

# Prune partial time directories written by the preCICE adapter (known issue #26).
echo "→ Pruning partial time directories (no H2O2 field)..."
for d in [0-9]* [0-9]*.[0-9]*; do
    [[ -d "$d" ]] || continue
    if [[ ! -f "$d/H2O2" ]]; then
        echo "  removed: $d"
        rm -rf "$d"
    fi
done

echo "→ Removing processor directories..."
rm -rf processor*/

touch tissue_42_decay.foam
echo "→ Done. Open tissue_42_decay.foam in ParaView."
