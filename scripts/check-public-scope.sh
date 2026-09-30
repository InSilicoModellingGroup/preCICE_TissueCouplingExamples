#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd "$(dirname "$0")/.." && pwd)
cd "$ROOT"

fail() {
  echo "scope check failed: $*" >&2
  exit 1
}

for directory in \
  adapters/openfoam-adapter/CHT \
  adapters/openfoam-adapter/FSI \
  adapters/openfoam-adapter/FF; do
  [[ ! -e "$directory" ]] || fail "excluded adapter module exists: $directory"
done

if rg -n 'UTrans|UCoeff|rhokTrans|ParticleMomentumSource' \
    adapters participants/cells/src; then
  fail "excluded particle feedback identifier found"
fi

if rg -n '/home/[^/]+/' . \
    --glob '!scripts/check-public-scope.sh' \
    --glob '!*.msh'; then
  fail "machine-specific home path found"
fi

expected_cases='tissue_41_diffusion tissue_42_decay tissue_43_cell_response tissue_44_twoway'
actual_cases=$(find participants -mindepth 1 -maxdepth 1 \
  -type d ! -name cells -printf '%f\n' | sort | tr '\n' ' ' | sed 's/ $//')
[[ "$actual_cases" == "$expected_cases" ]] || \
  fail "unexpected participant set: $actual_cases"

if find participants/tissue_4* -mindepth 1 -maxdepth 1 \
    -type d -regextype posix-extended -regex '.*/[1-9][0-9.]*' | grep -q .; then
  fail "generated OpenFOAM time directory found"
fi

if find . -type d \( -name processor0 -o -name postProcessing -o -name output \
    -o -name precice-profiling \) | grep -q .; then
  fail "generated simulation directory found"
fi

echo "public-scope check passed"
