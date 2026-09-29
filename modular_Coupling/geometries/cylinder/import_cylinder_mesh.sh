#!/bin/bash
# Import shared cylinder.msh into an OpenFOAM case.
# Accepts Gmsh 4.x; converts to a temporary msh2 file for gmshToFoam.
# The canonical mesh remains in this shared geometry directory; cases receive
# only the generated OpenFOAM polyMesh, which is ignored by Git.
# Usage: import_cylinder_mesh.sh [case_dir]
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
CASE_DIR=$(cd "${1:-$(pwd)}" && pwd)
MSH_SRC="${CYLINDER_MSH:-$SCRIPT_DIR/cylinder.msh}"

if [[ ! -f "$MSH_SRC" ]]; then
  echo "ERROR: cylinder mesh not found: $MSH_SRC" >&2
  exit 1
fi

set +u
source /usr/lib/openfoam/openfoam2406/etc/bashrc
set -euo pipefail

cd "$CASE_DIR"
mkdir -p constant system logs

# Prefer msh2 for OpenFOAM gmshToFoam, but never keep that conversion as
# another geometry source. This prevents stale meshes from surviving updates.
MSH2=$(mktemp "$CASE_DIR/logs/cylinder_msh2.XXXXXX.msh")
trap 'rm -f "$MSH2"' EXIT
echo "→ Converting $MSH_SRC to temporary msh2..."
gmsh -3 -format msh2 "$MSH_SRC" -o "$MSH2" > "$CASE_DIR/logs/gmsh_convert.log" 2>&1

rm -f constant/cylinder_msh2.msh constant/cylinder.vtk constant/cylinder.msh
rm -rf constant/polyMesh

echo "→ gmshToFoam temporary msh2 conversion"
gmshToFoam "$MSH2" > logs/gmshToFoam.log 2>&1
tail -30 logs/gmshToFoam.log

if [[ -f constant/polyMesh/boundary ]]; then
  sed -i -E 's/type[[:space:]]+patch;/type            wall;/' constant/polyMesh/boundary
  sed -i '/physicalType/d' constant/polyMesh/boundary
fi

echo "→ Patches:"
grep -E '^[ ]+[a-zA-Z]|nFaces' constant/polyMesh/boundary | head -40

echo "→ checkMesh"
checkMesh > logs/checkMesh.log 2>&1
grep -E 'bounding|cells:|Mesh OK|Failed' logs/checkMesh.log || true
tail -8 logs/checkMesh.log

echo "→ Cylinder mesh imported into $CASE_DIR"
