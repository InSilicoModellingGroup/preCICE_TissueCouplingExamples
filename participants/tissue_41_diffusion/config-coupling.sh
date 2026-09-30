#!/bin/bash

SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO_ROOT=$(cd "$SCRIPT_DIR/../.." && pwd)

set +u
source /usr/lib/openfoam/openfoam2406/etc/bashrc
set -euo pipefail
export LD_LIBRARY_PATH="$REPO_ROOT/adapters/openfoam-adapter/Make/$WM_OPTIONS:${LD_LIBRARY_PATH:-}"

cd "$SCRIPT_DIR"

sed "s|EXCHANGE_DIR|$REPO_ROOT/precice-run|g" \
    precice-config-tissue_41_diffusion.xml.in \
    > precice-config-tissue_41_diffusion.xml

echo "=========================================="
echo "  tissue_41_diffusion: Configure Coupled Run"
echo "=========================================="

foamDictionary -entry functions.preCICE_Adapter.enabled -set true  system/controlDict
foamDictionary -entry functions.runTimeControl.enabled  -set false system/controlDict
foamDictionary -entry purgeWrite -set 0 system/controlDict

if [ ! -f system/decomposeParDict ]; then
    echo "ERROR: Missing system/decomposeParDict" >&2
    exit 1
fi

NPROCS=$(awk '/numberOfSubdomains/{ gsub(";",""); print $2; exit }' system/decomposeParDict)
if [ -z "$NPROCS" ] || [ "$NPROCS" -le 1 ]; then
    echo "ERROR: tissue_41_diffusion coupling requires numberOfSubdomains > 1 in system/decomposeParDict" >&2
    exit 1
fi

if [[ ! -f constant/polyMesh/points ]]; then
    echo "→ Importing cylinder mesh..."
    ../../inputs/cylinder/import_cylinder_mesh.sh .
fi

if compgen -G "processor*" >/dev/null; then
    HAS_RESULT_TIME=false
    if [[ -d processor0 ]]; then
        for d in processor0/*; do
            [[ -d "$d" ]] || continue
            time_name=${d##*/}
            if [[ "$time_name" =~ ^[0-9]+([.][0-9]+)?$ ]] &&
               awk -v t="$time_name" 'BEGIN { exit !(t > 0) }'; then
                HAS_RESULT_TIME=true
                break
            fi
        done
    fi
    if [[ "$HAS_RESULT_TIME" == true ]]; then
        echo "→ Existing processor directories found; reconstructing before changing decomposition..."
        mkdir -p logs
        reconstructPar 2>&1 | tee logs/tissue_41_diffusion.config-reconstruct.log
    else
        echo "→ Existing processor directories contain no result times; skipping reconstruction..."
    fi
    echo "→ Removing old processor directories..."
    rm -rf processor*/
fi

echo "→ Decomposing with $NPROCS ranks..."
decomposePar -force

ACTUAL_NPROCS=$(find processor* -maxdepth 0 -type d 2>/dev/null | wc -l)
if [ "$ACTUAL_NPROCS" -ne "$NPROCS" ]; then
    echo "ERROR: Expected $NPROCS processor dirs, found $ACTUAL_NPROCS." >&2
    exit 1
fi

if [[ ! -f processor0/0/H2O2 || ! -f processor0/0/NO2 ]]; then
    echo "ERROR: decomposePar did not write processor*/0/{H2O2,NO2}." >&2
    echo "       Check that 0/H2O2 and 0/NO2 have valid FoamFile headers." >&2
    exit 1
fi

echo ""
echo "=========================================="
echo "  tissue_41_diffusion Coupled Run Configured"
echo "=========================================="
echo "  Ranks: $NPROCS"
echo ""
echo "Next step:"
echo "  ./run.sh"
