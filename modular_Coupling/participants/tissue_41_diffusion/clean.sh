#!/bin/bash

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
RUN_LOCK_DIR="$SCRIPT_DIR/.tissue_41_diffusion-run.lock"

list_case_solver_pids() {
	local pid
	local cwd
	local cmdline

	if ! command -v pgrep >/dev/null 2>&1; then
		return
	fi

	while IFS= read -r pid; do
		[[ "$pid" =~ ^[0-9]+$ ]] || continue
		cwd="$(readlink -f "/proc/$pid/cwd" 2>/dev/null || true)"
		[[ "$cwd" == "$SCRIPT_DIR" ]] || continue

		cmdline="$(tr '\0' ' ' < "/proc/$pid/cmdline" 2>/dev/null || true)"
		if [[ "$cmdline" == *"myDiffusionFoam"* || "$cmdline" == *"mpirun"* ]]; then
			printf '%s\n' "$pid"
		fi
	done < <(pgrep -f 'myDiffusionFoam|mpirun' 2>/dev/null || true)
}

print_case_solver_processes() {
	local pid
	local cmdline
	for pid in "$@"; do
		cmdline="$(tr '\0' ' ' < "/proc/$pid/cmdline" 2>/dev/null || true)"
		echo "         pid=$pid cmd=${cmdline:-<unavailable>}"
	done
}

ensure_no_active_case_solver() {
	local -a active_pids
	mapfile -t active_pids < <(list_case_solver_pids | sort -n -u)

	if [[ "${#active_pids[@]}" -eq 0 ]]; then
		return
	fi

	echo "ERROR: Active OpenFOAM process(es) detected for this case directory." >&2
	print_case_solver_processes "${active_pids[@]}" >&2
	echo "       Stop these processes before cleaning." >&2
	exit 1
}

cd "$SCRIPT_DIR"

if [[ "${1:-}" == "--help" ]]; then
	echo "Usage: ./clean.sh"
	echo "Removes generated tissue_41_diffusion OpenFOAM outputs and resets controlDict."
	exit 0
fi

if [[ "$#" -gt 0 ]]; then
	echo "ERROR: Unknown option: $1" >&2
	echo "Run ./clean.sh --help for usage." >&2
	exit 1
fi

echo "Cleaning OpenFOAM simulation output in: $SCRIPT_DIR"

ensure_no_active_case_solver

if [[ -f "$RUN_LOCK_DIR/pid" ]]; then
	run_pid="$(tr -d '[:space:]' < "$RUN_LOCK_DIR/pid" 2>/dev/null || true)"
	if [[ "$run_pid" =~ ^[0-9]+$ ]] && kill -0 "$run_pid" 2>/dev/null; then
		echo "ERROR: Active tissue_41_diffusion run detected (pid=$run_pid)." >&2
		exit 1
	fi
	echo "WARNING: Removing stale tissue_41_diffusion run lock at $RUN_LOCK_DIR"
fi

echo "Removing time directories..."
find . -mindepth 1 -maxdepth 1 -type d -printf '%f\n' \
	| LC_ALL=C awk '$1 ~ /^[0-9]+(\.[0-9]+)?$/ && ($1 + 0) > 0 {print $1}' \
	| xargs -r rm -rf

echo "Resetting controlDict to initial state..."
if [[ -f system/controlDict ]]; then
	sed -i -E 's/^[[:space:]]*startTime[[:space:]]+[0-9.eE+-]+;/startTime       0;/' system/controlDict
	sed -i -E 's/^[[:space:]]*endTime[[:space:]]+[0-9.eE+-]+;/endTime         180;/' system/controlDict
	sed -i -E 's/^[[:space:]]*startFrom[[:space:]]+latestTime;/startFrom       startTime;/' system/controlDict
fi

echo "Preserving validated cylinder mesh (constant/polyMesh)."

echo "Removing log files..."
rm -f tissue_41_diffusion.log
rm -f log.*
rm -rf logs/*
rm -f precice-*-iterations.log
rm -f precice-*-convergence.log
rm -f precice-*-watchpoint-*.log
rm -f precice-*-watchintegral-*.log
rm -f profiling.json trace.json core
rm -rf precice-exports

if [[ -d precice-profiling ]]; then
	rm -rf precice-profiling/*
fi

echo "Removing functionObject and post-processing artifacts..."
rm -rf postProcessing history dynamicCode
rm -rf 0/uniform/functionObjects/functionObjectProperties

echo "Removing processor directories..."
rm -rf processor*/

rm -rf "$RUN_LOCK_DIR"
rm -f p*_times.txt all_possible_times.txt

echo "Removing other generated files..."
rm -f *.dat *.vtk *.foam *.csv *~

echo "Clean up complete."
