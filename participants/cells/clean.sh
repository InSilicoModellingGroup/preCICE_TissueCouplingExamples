#!/bin/bash

# Script to clean up BioDynaMo generated files in the cells directory

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO_ROOT=$(cd "$SCRIPT_DIR/../.." && pwd)
PRECICE_RUN_DIR="$REPO_ROOT/precice-run"

cd "$SCRIPT_DIR"

echo "Cleaning BioDynaMo simulation output..."

# Remove log files
echo "Removing log files..."
rm -f cells.log
rm -f precice-cells-convergence.log precice-cells-iterations.log

# Clean output directory but keep the directory itself
if [ -d "output" ]; then
  echo "Cleaning output directory..."
  rm -rf output/*
fi

# Remove precice-profiling directory contents
if [ -d "precice-profiling" ]; then
  echo "Cleaning precice-profiling directory..."
  rm -rf precice-profiling/
fi
rm -rf "$PRECICE_RUN_DIR"

# Remove symlinked precice-config.xml (created by run.sh)
if [ -L "precice-config.xml" ]; then
  echo "Removing precice-config.xml symlink..."
  rm -f precice-config.xml
fi
if [ -L "params.csv" ]; then
  echo "Removing params.csv symlink..."
  rm -f params.csv
fi

# Clean build artifacts completely
if [ -d "build" ]; then
  echo "Removing build directory..."
  rm -rf build/
fi

# Remove any temporary files
echo "Removing temporary files..."
rm -f *~
rm -f *.dat
rm -f core.*

echo "Clean up complete."
