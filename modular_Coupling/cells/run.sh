#!/bin/bash

# set -euo pipefail

# Simple wrapper to source BioDynaMo and run the simulation
# Usage: ./run.sh [participant_name]
# Example: ./run.sh simple
#          ./run.sh diff
#          ./run.sh conv_diff
#          ./run.sh venturi

# Determine script directory
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
MODULAR_COUPLING_DIR=$(cd "$SCRIPT_DIR/.." && pwd)
export HWLOC_HIDE_ERRORS="${HWLOC_HIDE_ERRORS:-2}"

# Parse participant name from first argument (default: simple)
PARTICIPANT="${1:-simple}"
shift 1 2>/dev/null || true

# Validate participant exists
PARTICIPANT_DIR="$MODULAR_COUPLING_DIR/participants/$PARTICIPANT"
if [[ ! -d "$PARTICIPANT_DIR" ]]; then
  echo "Error: Participant '$PARTICIPANT' not found in participants/" >&2
  echo "Available participants:"
  ls -1 "$MODULAR_COUPLING_DIR/participants/" 2>/dev/null | sed 's/^/  - /'
  exit 1
fi

echo "=========================================="
echo "Coupling with participant: $PARTICIPANT"
echo "=========================================="

if [[ "$PARTICIPANT" == "o2Cube" ]]; then
  echo "Role info: cells is the second participant/socket connector for o2Cube."
  echo "Start OpenFOAM first: cd ../participants/o2Cube && ./run.sh"
fi

if [[ "$PARTICIPANT" == "tissue_modelling" || "$PARTICIPANT" == tissue_4* ]]; then
  echo "Role info: cells is the second participant/socket connector for $PARTICIPANT."
  echo "Start OpenFOAM first: cd ../participants/$PARTICIPANT && ./run.sh"
fi

# Select the correct params CSV for this participant
PARAMS_CONFIG="$SCRIPT_DIR/params-$PARTICIPANT.csv"
if [[ -f "$PARAMS_CONFIG" ]]; then
  echo "Using params config: $PARAMS_CONFIG"
else
  PARAMS_CONFIG="$SCRIPT_DIR/params.csv"
  if [[ -f "$PARAMS_CONFIG" ]]; then
    echo "Using default params config: $PARAMS_CONFIG"
  else
    echo "Error: params file not found for participant '$PARTICIPANT'" >&2
    exit 1
  fi
fi

RUN_PARAMS_CONFIG="$PARAMS_CONFIG"
echo "Run params: $RUN_PARAMS_CONFIG"

# Generate the runtime preCICE config from the tracked, portable template.
PRECICE_CONFIG_REL=$(awk -F',' 'tolower($1) ~ /precice\.config_file/ {gsub(/^[ \t]+|[ \t]+$/, "", $3); print $3; exit}' "$RUN_PARAMS_CONFIG")
if [[ -n "$PRECICE_CONFIG_REL" ]]; then
  if [[ "$PRECICE_CONFIG_REL" = /* ]]; then
    PRECICE_CONFIG_PATH="$PRECICE_CONFIG_REL"
  else
    PRECICE_CONFIG_PATH="$SCRIPT_DIR/$PRECICE_CONFIG_REL"
  fi
  PRECICE_CONFIG_TEMPLATE="$PRECICE_CONFIG_PATH.in"
  if [[ -f "$PRECICE_CONFIG_TEMPLATE" ]]; then
    EXCHANGE_DIR="$MODULAR_COUPLING_DIR/precice-run"
    mkdir -p "$EXCHANGE_DIR"
    sed "s|EXCHANGE_DIR|$EXCHANGE_DIR|g" \
      "$PRECICE_CONFIG_TEMPLATE" > "$PRECICE_CONFIG_PATH"
    echo "Generated runtime preCICE config: $PRECICE_CONFIG_PATH"
  fi
fi

# This report repository does not vendor BioDynaMo. Point BDM_PATH to its
# thisbdm.sh unless an active BioDynaMo environment is already available.
BDM_PATH="${BDM_PATH:-}"

# Ensure BioDynaMo helper variables exist even with nounset enabled
: "${BDM_THISBDM_LOGLEVEL:=5}"
: "${BDM_THISBDM_QUIET:=false}"
: "${BDM_THISBDM_SILENT:=false}"

# Source BioDynaMo environment
if [[ -z "${BDMSYS:-}" ]]; then
  if [[ -z "$BDM_PATH" || ! -f "$BDM_PATH" ]]; then
    echo "Error: set BDM_PATH to BioDynaMo's thisbdm.sh." >&2
    exit 1
  fi
  echo "Sourcing BioDynaMo environment from $BDM_PATH"
  source "$BDM_PATH"
fi

# Check if BioDynaMo is available
if [[ -z "${BDMSYS:-}" ]]; then
  echo "Error: Could not load BioDynaMo environment." >&2
  exit 1
fi

echo "Using BioDynaMo from: $BDMSYS"

# Build and run the simulation directly.
echo "Building and running BioDynaMo simulation..."
cd "$SCRIPT_DIR"

# Set OpenFOAM case directory based on participant
OF_CASE_DIR="$PARTICIPANT_DIR"
STEADY_FILE="$OF_CASE_DIR/.steady_state_time"
CONTROL_DICT="$OF_CASE_DIR/system/controlDict"

# Export participant name for BioDynaMo to use
export BDM_PARTICIPANT="$PARTICIPANT"

read_params_entry() {
  local key="$1"
  if [[ ! -f "$RUN_PARAMS_CONFIG" ]]; then
    return
  fi
  awk -F',' -v entry="$key" '
    NR == 1 { next }
    {
      name = $1
      value = $3
      gsub(/^[ \t]+|[ \t]+$/, "", name)
      gsub(/^[ \t]+|[ \t]+$/, "", value)
      if (name == entry) {
        print value
        exit
      }
    }
  ' "$RUN_PARAMS_CONFIG" 2>/dev/null
}

resolve_case_output_dir() {
  local configured_dir
  configured_dir="$(read_params_entry "visualization.output_directory")"

  if [[ -z "$configured_dir" ]]; then
    configured_dir="output/$PARTICIPANT"
  fi

  if [[ "$configured_dir" = /* ]]; then
    printf '%s\n' "$configured_dir"
  else
    printf '%s\n' "$SCRIPT_DIR/$configured_dir"
  fi
}

collect_precice_artifacts() {
  local src

  for src in "precice-cells-convergence.log" "precice-cells-iterations.log"; do
    if [[ -f "$SCRIPT_DIR/$src" ]]; then
      mv -f "$SCRIPT_DIR/$src" "$CASE_LOG_DIR/$src"
    fi
  done

  if [[ -d "$SCRIPT_DIR/precice-profiling" ]]; then
    mkdir -p "$CASE_PRECICE_PROFILE_DIR"
    if compgen -G "$SCRIPT_DIR/precice-profiling/*" >/dev/null; then
      mv -f "$SCRIPT_DIR/precice-profiling"/* "$CASE_PRECICE_PROFILE_DIR/" 2>/dev/null || true
    fi
    rmdir "$SCRIPT_DIR/precice-profiling" 2>/dev/null || true
  fi
}

# Ensure preCICE profiling is off when disabled in params (avoid stale env).
PARAM_PRECICE_PROFILING="$(read_params_entry "precice.enable_profiling")"
if [[ "${PARAM_PRECICE_PROFILING,,}" == "false" || "$PARAM_PRECICE_PROFILING" == "0" ]]; then
  unset PRECICE_PROFILING
  unset PRECICE_PROFILING_MODE
fi

CASE_OUTPUT_DIR="$(resolve_case_output_dir)"
CASE_LOG_DIR="$CASE_OUTPUT_DIR/logs"
CASE_PRECICE_PROFILE_DIR="$CASE_OUTPUT_DIR/precice-profiling"
LOG_FILE="$CASE_LOG_DIR/cells_${PARTICIPANT}.log"

mkdir -p "$CASE_OUTPUT_DIR" "$CASE_LOG_DIR" "$CASE_PRECICE_PROFILE_DIR"

echo "Case output directory: $CASE_OUTPUT_DIR"
echo "Case log directory:    $CASE_LOG_DIR"

# The bundled BioDynaMo PVTU writer names files by scheduler step.  Pass the
# physical seconds per scheduler step to the ParaView-state generator so it can
# build a PVD collection with an actual 0..max-time clock.
BDM_PHYSICAL_TIME_STEP="$(read_params_entry "simulation.biodynamo_timestep")"
if [[ -n "$BDM_PHYSICAL_TIME_STEP" ]]; then
  export BDM_PHYSICAL_TIME_STEP
  echo "BioDynaMo visualization time scale: $BDM_PHYSICAL_TIME_STEP s/step"
fi

# Reset participant-local logs/profiling snapshots for a fresh run.
rm -f "$LOG_FILE"
rm -f "$CASE_LOG_DIR/precice-cells-convergence.log" "$CASE_LOG_DIR/precice-cells-iterations.log"
rm -rf "$CASE_PRECICE_PROFILE_DIR"/*

# Remove root-level preCICE logs so they can be re-collected into the case directory.
rm -f "$SCRIPT_DIR/precice-cells-convergence.log" "$SCRIPT_DIR/precice-cells-iterations.log"
if [[ -d "$SCRIPT_DIR/precice-profiling" ]]; then
  rm -rf "$SCRIPT_DIR/precice-profiling"/*
fi

detect_total_logical_cpus() {
  if command -v nproc >/dev/null 2>&1; then
    nproc
    return
  fi
  getconf _NPROCESSORS_ONLN 2>/dev/null || echo 1
}

detect_physical_cpu_list() {
  if command -v lscpu >/dev/null 2>&1; then
    local parsed
    parsed="$({ lscpu -p=CPU,CORE,SOCKET 2>/dev/null || true; } | awk -F',' '
      /^#/ { next }
      {
        cpu = $1 + 0
        core = $2
        socket = $3
        key = socket ":" core
        if (!(key in min_cpu) || cpu < min_cpu[key]) {
          min_cpu[key] = cpu
        }
      }
      END {
        for (k in min_cpu) {
          print min_cpu[k]
        }
      }
    ' | sort -n)"
    if [[ -n "$parsed" ]]; then
      printf '%s\n' "$parsed"
      return
    fi
  fi

  local logical
  logical="$(detect_total_logical_cpus)"
  local physical=$((logical / 2))
  if [[ "$physical" -lt 1 ]]; then
    physical=1
  fi
  seq 0 $((physical - 1))
}

detect_logical_cpu_list() {
  if command -v lscpu >/dev/null 2>&1; then
    local parsed
    parsed="$({ lscpu -p=CPU 2>/dev/null || true; } | awk -F',' '
      /^#/ { next }
      { print $1 + 0 }
    ' | sort -n)"
    if [[ -n "$parsed" ]]; then
      printf '%s\n' "$parsed"
      return
    fi
  fi

  local logical
  logical="$(detect_total_logical_cpus)"
  seq 0 $((logical - 1))
}

detect_sibling_cpuset_for_cpu() {
  local target_cpu="${1:-0}"

  if command -v lscpu >/dev/null 2>&1; then
    local parsed
    parsed="$({ lscpu -p=CPU,CORE,SOCKET 2>/dev/null || true; } | awk -F',' -v target="$target_cpu" '
      /^#/ { next }
      {
        cpu = $1 + 0
        core = $2
        socket = $3
        if (cpu == target) {
          target_core = core
          target_socket = socket
        }
        cpu_list[cpu] = core "," socket
      }
      END {
        for (cpu in cpu_list) {
          split(cpu_list[cpu], parts, ",")
          if (parts[1] == target_core && parts[2] == target_socket) {
            print cpu
          }
        }
      }
    ' | sort -n | paste -sd, -)"
    if [[ -n "$parsed" ]]; then
      printf '%s\n' "$parsed"
      return
    fi
  fi

  printf '%s\n' "$target_cpu"
}

detect_physical_anchor_cpu_for_cpu() {
  local target_cpu="${1:-0}"

  if command -v lscpu >/dev/null 2>&1; then
    local parsed
    parsed="$({ lscpu -p=CPU,CORE,SOCKET 2>/dev/null || true; } | awk -F',' -v target="$target_cpu" '
      /^#/ { next }
      {
        cpu = $1 + 0
        core = $2
        socket = $3
        key = socket ":" core
        if (!(key in min_cpu) || cpu < min_cpu[key]) {
          min_cpu[key] = cpu
        }
        if (cpu == target) {
          target_key = key
        }
      }
      END {
        if (target_key in min_cpu) {
          print min_cpu[target_key]
        }
      }
    ' | tail -n1)"
    if [[ -n "$parsed" ]]; then
      printf '%s\n' "$parsed"
      return
    fi
  fi

  printf '%s\n' "$target_cpu"
}

build_cpuset_excluding_reserved() {
  local reserved_cpuset="$1"
  local requested_count="$2"

  detect_logical_cpu_list | awk -v reserved="$reserved_cpuset" -v limit="$requested_count" '
    BEGIN {
      reserved_count = split(reserved, tokens, ",")
      for (i = 1; i <= reserved_count; ++i) {
        token = tokens[i]
        gsub(/[[:space:]]/, "", token)
        if (token == "") {
          continue
        }
        if (token ~ /-/) {
          split(token, bounds, "-")
          start = bounds[1] + 0
          stop = bounds[2] + 0
          for (cpu = start; cpu <= stop; ++cpu) {
            excluded[cpu] = 1
          }
        } else {
          excluded[token + 0] = 1
        }
      }
      selected = 0
    }
    {
      cpu = $1 + 0
      if (!(cpu in excluded)) {
        cpus[selected] = cpu
        selected += 1
        if (selected >= limit) {
          exit
        }
      }
    }
    END {
      for (i = 0; i < selected; ++i) {
        printf "%s%d", (i > 0 ? "," : ""), cpus[i]
      }
      printf "\n"
    }
  '
}

count_cpuset_entries() {
  local cpuset="$1"
  awk -F',' '
    {
      count = 0
      for (i = 1; i <= NF; ++i) {
        token = $i
        gsub(/[[:space:]]/, "", token)
        if (token == "") {
          continue
        }
        if (token ~ /-/) {
          split(token, bounds, "-")
          start = bounds[1] + 0
          stop = bounds[2] + 0
          if (stop >= start) {
            count += (stop - start + 1)
          }
        } else {
          count += 1
        }
      }
      print count
    }
  ' <<< "$cpuset"
}

build_cpuset_from_slice() {
  local start="$1"
  local length="$2"
  local -a slice=("${PHYSICAL_CPUS[@]:start:length}")
  local IFS=','
  echo "${slice[*]}"
}

build_cpuset_excluding_physical_cpu() {
  local excluded_cpu="$1"
  local length="$2"
  local anchor_cpu
  local -a selected=()
  local cpu
  local IFS=','

  anchor_cpu="$(detect_physical_anchor_cpu_for_cpu "$excluded_cpu")"

  for cpu in "${PHYSICAL_CPUS[@]}"; do
    if [[ "$cpu" == "$anchor_cpu" ]]; then
      continue
    fi
    selected+=("$cpu")
    if [[ "${#selected[@]}" -ge "$length" ]]; then
      break
    fi
  done

  echo "${selected[*]}"
}

build_physical_cpuset_excluding_reserved() {
  local reserved_cpuset="$1"
  local requested_count="$2"
  local physical_csv

  local IFS=','
  physical_csv="${PHYSICAL_CPUS[*]}"

  awk -v cpus="$physical_csv" -v reserved="$reserved_cpuset" -v limit="$requested_count" '
    BEGIN {
      cpu_count = split(cpus, cpu_tokens, ",")

      reserved_count = split(reserved, tokens, ",")
      for (i = 1; i <= reserved_count; ++i) {
        token = tokens[i]
        gsub(/[[:space:]]/, "", token)
        if (token == "") {
          continue
        }
        if (token ~ /-/) {
          split(token, bounds, "-")
          start = bounds[1] + 0
          stop = bounds[2] + 0
          for (cpu = start; cpu <= stop; ++cpu) {
            excluded[cpu] = 1
          }
        } else {
          excluded[token + 0] = 1
        }
      }

      selected = 0
      for (i = 1; i <= cpu_count; ++i) {
        cpu = cpu_tokens[i] + 0
        if (!(cpu in excluded)) {
          chosen[selected] = cpu
          selected += 1
          if (selected >= limit) {
            break
          }
        }
      }

      for (i = 0; i < selected; ++i) {
        printf "%s%d", (i > 0 ? "," : ""), chosen[i]
      }
      printf "\n"
    }
  '
}

read_participant_subdomains() {
  local dict="$PARTICIPANT_DIR/system/decomposeParDict"
  if [[ -f "$dict" ]]; then
    awk '/numberOfSubdomains/{gsub(";","",$2); print $2; exit}' "$dict" 2>/dev/null
  fi
}

detect_openfoam_layout_mode() {
  case "$PARTICIPANT" in
    simple|diff|venturi_a)
      printf '%s\n' "serial"
      ;;
    conv_diff|o2Cube|tissue_modelling|tissue_41_diffusion|tissue_42_decay|tissue_43_cell_response|tissue_44_twoway|aneurysm|venturi_b|venturi_c)
      printf '%s\n' "parallel"
      ;;
    *)
      printf '%s\n' "unknown"
      ;;
  esac
}

detect_openfoam_parallel_ranks() {
  case "$PARTICIPANT" in
    venturi_b)
      printf '%s\n' "${VENTURI_B_COUPLING_NPROCS:-10}"
      ;;
    venturi_c)
      printf '%s\n' "${VENTURI_C_COUPLING_NPROCS:-10}"
      ;;
    conv_diff|o2Cube|tissue_modelling|tissue_41_diffusion|tissue_42_decay|tissue_43_cell_response|tissue_44_twoway|aneurysm)
      read_participant_subdomains
      ;;
    *)
      printf '%s\n' "1"
      ;;
  esac
}

detect_openfoam_serial_pin_cpu() {
  case "$PARTICIPANT" in
    simple)
      printf '%s\n' "${SIMPLE_OF_PIN_CPU:-0}"
      ;;
    diff)
      printf '%s\n' "${DIFF_OF_PIN_CPU:-0}"
      ;;
    venturi_a)
      printf '%s\n' "${VENTURI_A_OF_PIN_CPU:-0}"
      ;;
    *)
      printf '%s\n' "0"
      ;;
  esac
}

# Configure OpenMP parallelization for BioDynaMo.
# Thread count is centrally managed here from physical-core availability
# after accounting for OpenFOAM MPI ranks and OS reserve.
OF_SUBDOMAINS=""
RESERVED_SYSTEM_CPUS="${BDM_RESERVED_SYSTEM_CPUS:-5}"
if [[ ! "$RESERVED_SYSTEM_CPUS" =~ ^[0-9]+$ ]]; then
  RESERVED_SYSTEM_CPUS=5
fi

mapfile -t PHYSICAL_CPUS < <(detect_physical_cpu_list)
PHYSICAL_CORE_COUNT="${#PHYSICAL_CPUS[@]}"
if [[ "$PHYSICAL_CORE_COUNT" -lt 1 ]]; then
  PHYSICAL_CORE_COUNT=1
  PHYSICAL_CPUS=(0)
fi

OF_LAYOUT_MODE="$(detect_openfoam_layout_mode)"
OF_PHYSICAL_CORES=0
OF_SERIAL_PIN_CPU=""
if [[ "$OF_LAYOUT_MODE" == "serial" ]]; then
  OF_PHYSICAL_CORES=1
  OF_SERIAL_PIN_CPU="$(detect_openfoam_serial_pin_cpu)"
elif [[ "$OF_LAYOUT_MODE" == "parallel" ]]; then
  OF_PHYSICAL_CORES="$(detect_openfoam_parallel_ranks)"
  if [[ ! "$OF_PHYSICAL_CORES" =~ ^[0-9]+$ || "$OF_PHYSICAL_CORES" -lt 1 ]]; then
    OF_PHYSICAL_CORES=1
  fi
fi

SAFE_AVAILABLE_PHYSICAL=$((PHYSICAL_CORE_COUNT - OF_PHYSICAL_CORES - RESERVED_SYSTEM_CPUS))
if [[ "$SAFE_AVAILABLE_PHYSICAL" -lt 1 ]]; then
  SAFE_AVAILABLE_PHYSICAL=1
fi

if [[ -n "${BDM_OMP_THREADS:-}" ]]; then
  OMP_THREADS_SELECTED="$BDM_OMP_THREADS"
  OMP_THREADS_SOURCE="env(BDM_OMP_THREADS)"
elif [[ "$OF_LAYOUT_MODE" == "serial" || "$OF_LAYOUT_MODE" == "parallel" ]]; then
  OMP_THREADS_SELECTED="$SAFE_AVAILABLE_PHYSICAL"
  OMP_THREADS_SOURCE="auto(physical_cores=$PHYSICAL_CORE_COUNT - openfoam_cores=$OF_PHYSICAL_CORES - reserved_system_cores=$RESERVED_SYSTEM_CPUS)"
else
  OMP_THREADS_SELECTED=30
  OMP_THREADS_SOURCE="default(30)"
fi

if [[ ! "$OMP_THREADS_SELECTED" =~ ^[0-9]+$ || "$OMP_THREADS_SELECTED" -lt 1 ]]; then
  OMP_THREADS_SELECTED=1
  OMP_THREADS_SOURCE="$OMP_THREADS_SOURCE, corrected_to_minimum(1)"
fi

if [[ -n "${BDM_OMP_MAX_THREADS:-}" && "${BDM_OMP_MAX_THREADS}" =~ ^[0-9]+$ && "${BDM_OMP_MAX_THREADS}" -gt 0 ]]; then
  if [[ "$OMP_THREADS_SELECTED" -gt "$BDM_OMP_MAX_THREADS" ]]; then
    OMP_THREADS_SELECTED="$BDM_OMP_MAX_THREADS"
    OMP_THREADS_SOURCE="$OMP_THREADS_SOURCE, capped_by_env(BDM_OMP_MAX_THREADS=$BDM_OMP_MAX_THREADS)"
  fi
fi

export OMP_NUM_THREADS="$OMP_THREADS_SELECTED"
export OMP_DYNAMIC=false
export OMP_PROC_BIND=spread
export OMP_PLACES=cores
export OMP_NESTED=false
export OMP_MAX_ACTIVE_LEVELS=1
echo "BioDynaMo OpenMP configuration: OMP_NUM_THREADS=$OMP_NUM_THREADS (source=$OMP_THREADS_SOURCE) OMP_DYNAMIC=$OMP_DYNAMIC OMP_PROC_BIND=$OMP_PROC_BIND OMP_PLACES=$OMP_PLACES"

BDM_CPUSET="${BDM_CPUSET:-}"
if [[ "$OF_LAYOUT_MODE" == "serial" ]]; then
  OF_RESERVED_CPUSET="$(detect_sibling_cpuset_for_cpu "$OF_SERIAL_PIN_CPU")"
  if [[ -z "$BDM_CPUSET" ]]; then
    BDM_CPUSET="$(build_cpuset_excluding_physical_cpu "$OF_SERIAL_PIN_CPU" "$SAFE_AVAILABLE_PHYSICAL")"
  fi

  bdm_cpuset_count="$(count_cpuset_entries "$BDM_CPUSET")"
  if [[ "$bdm_cpuset_count" =~ ^[0-9]+$ && "$bdm_cpuset_count" -gt 0 && "$OMP_NUM_THREADS" -gt "$bdm_cpuset_count" ]]; then
    export OMP_NUM_THREADS="$bdm_cpuset_count"
    echo "Adjusted OMP_NUM_THREADS to match BDM_CPUSET capacity: OMP_NUM_THREADS=$OMP_NUM_THREADS"
  elif [[ ! "$bdm_cpuset_count" =~ ^[0-9]+$ || "$bdm_cpuset_count" -lt 1 ]]; then
    export OMP_NUM_THREADS=1
    echo "WARNING: No isolated BioDynaMo cpuset could be derived for $PARTICIPANT; falling back to OMP_NUM_THREADS=1 without explicit CPU pinning."
  fi

  export BDM_CPUSET

  echo "$PARTICIPANT physical-core layout:"
  echo "  physical_cores=$PHYSICAL_CORE_COUNT reserved_cores=$RESERVED_SYSTEM_CPUS"
  echo "  OpenFOAM pinned CPU: $OF_SERIAL_PIN_CPU"
  echo "  Reserved OpenFOAM sibling set: $OF_RESERVED_CPUSET"
  echo "  BioDynaMo cores (BDM_CPUSET): $BDM_CPUSET"
fi

if [[ "$OF_LAYOUT_MODE" == "parallel" ]]; then
  OF_CPUSET=""
  case "$PARTICIPANT" in
    o2Cube)
      OF_CPUSET="${O2CUBE_CPUSET:-$(build_cpuset_from_slice 0 "$OF_PHYSICAL_CORES")}"
      export O2CUBE_CPUSET="$OF_CPUSET"
      ;;
    tissue_modelling)
      OF_CPUSET="${TISSUE_MODELLING_CPUSET:-$(build_cpuset_from_slice 0 "$OF_PHYSICAL_CORES")}"
      export TISSUE_MODELLING_CPUSET="$OF_CPUSET"
      ;;
    tissue_41_diffusion)
      OF_CPUSET="${TISSUE_41_DIFFUSION_CPUSET:-$(build_cpuset_from_slice 0 "$OF_PHYSICAL_CORES")}"
      export TISSUE_41_DIFFUSION_CPUSET="$OF_CPUSET"
      ;;
    tissue_42_decay)
      OF_CPUSET="${TISSUE_42_DECAY_CPUSET:-$(build_cpuset_from_slice 0 "$OF_PHYSICAL_CORES")}"
      export TISSUE_42_DECAY_CPUSET="$OF_CPUSET"
      ;;
    tissue_43_cell_response)
      OF_CPUSET="${TISSUE_43_CELL_RESPONSE_CPUSET:-$(build_cpuset_from_slice 0 "$OF_PHYSICAL_CORES")}"
      export TISSUE_43_CELL_RESPONSE_CPUSET="$OF_CPUSET"
      ;;
    tissue_44_twoway)
      OF_CPUSET="${TISSUE_44_TWOWAY_CPUSET:-$(build_cpuset_from_slice 0 "$OF_PHYSICAL_CORES")}"
      export TISSUE_44_TWOWAY_CPUSET="$OF_CPUSET"
      ;;
    conv_diff)
      OF_CPUSET="${CONV_DIFF_OF_CPUSET:-$(build_cpuset_from_slice 0 "$OF_PHYSICAL_CORES")}"
      export CONV_DIFF_OF_CPUSET="$OF_CPUSET"
      ;;
    aneurysm)
      OF_CPUSET="${ANEURYSM_OF_CPUSET:-$(build_cpuset_from_slice 0 "$OF_PHYSICAL_CORES")}"
      export ANEURYSM_OF_CPUSET="$OF_CPUSET"
      ;;
    venturi_b)
      OF_CPUSET="${VENTURI_B_OF_CPUSET:-$(build_cpuset_from_slice 0 "$OF_PHYSICAL_CORES")}"
      export VENTURI_B_OF_CPUSET="$OF_CPUSET"
      ;;
    venturi_c)
      OF_CPUSET="${VENTURI_C_OF_CPUSET:-$(build_cpuset_from_slice 0 "$OF_PHYSICAL_CORES")}"
      export VENTURI_C_OF_CPUSET="$OF_CPUSET"
      ;;
  esac

  if [[ -z "$BDM_CPUSET" ]]; then
    if [[ -n "$OF_CPUSET" ]]; then
      BDM_CPUSET="$(build_physical_cpuset_excluding_reserved "$OF_CPUSET" "$SAFE_AVAILABLE_PHYSICAL")"
    else
      BDM_CPUSET="$(build_cpuset_from_slice "$OF_PHYSICAL_CORES" "$SAFE_AVAILABLE_PHYSICAL")"
    fi
  fi

  bdm_cpuset_count="$(count_cpuset_entries "$BDM_CPUSET")"
  if [[ "$bdm_cpuset_count" =~ ^[0-9]+$ && "$bdm_cpuset_count" -gt 0 && "$OMP_NUM_THREADS" -gt "$bdm_cpuset_count" ]]; then
    export OMP_NUM_THREADS="$bdm_cpuset_count"
    echo "Adjusted OMP_NUM_THREADS to match BDM_CPUSET capacity: OMP_NUM_THREADS=$OMP_NUM_THREADS"
  elif [[ ! "$bdm_cpuset_count" =~ ^[0-9]+$ || "$bdm_cpuset_count" -lt 1 ]]; then
    export OMP_NUM_THREADS=1
    echo "WARNING: No isolated BioDynaMo cpuset could be derived for $PARTICIPANT; falling back to OMP_NUM_THREADS=1 without explicit CPU pinning."
  fi
  export BDM_CPUSET

  RESERVED_CPUSET=""
  if [[ "$RESERVED_SYSTEM_CPUS" -gt 0 && -n "$OF_CPUSET" ]]; then
    RESERVED_CPUSET="$(build_physical_cpuset_excluding_reserved "$OF_CPUSET,$BDM_CPUSET" "$RESERVED_SYSTEM_CPUS")"
  fi

  echo "$PARTICIPANT physical-core layout:"
  echo "  physical_cores=$PHYSICAL_CORE_COUNT reserved_cores=$RESERVED_SYSTEM_CPUS"
  echo "  OpenFOAM physical cores: $OF_PHYSICAL_CORES"
  if [[ -n "$OF_CPUSET" ]]; then
    echo "  OpenFOAM cores: $OF_CPUSET"
  fi
  echo "  BioDynaMo cores (BDM_CPUSET): $BDM_CPUSET"
  if [[ -n "$RESERVED_CPUSET" ]]; then
    echo "  Reserved idle cores: $RESERVED_CPUSET"
  fi
  if [[ "$OF_PHYSICAL_CORES" -gt "$PHYSICAL_CORE_COUNT" ]]; then
    echo "WARNING: OpenFOAM requests $OF_PHYSICAL_CORES ranks but only $PHYSICAL_CORE_COUNT physical cores are available; strict conflict-free pinning is impossible without reducing OpenFOAM ranks."
  fi
fi

# Prevent hidden nested threading from BLAS/math libs during coupled runs.
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-1}"
echo "Aux thread caps: OPENBLAS_NUM_THREADS=$OPENBLAS_NUM_THREADS MKL_NUM_THREADS=$MKL_NUM_THREADS NUMEXPR_NUM_THREADS=$NUMEXPR_NUM_THREADS"

read_control_dict_entry() {
  local key="$1"
  if [[ ! -f "$CONTROL_DICT" ]]; then
    return
  fi
  awk -v entry="$key" '
    $1 == entry {
      gsub(";", "", $2);
      print $2;
      exit
    }
  ' "$CONTROL_DICT" 2>/dev/null
}

if [[ -f "$STEADY_FILE" ]]; then
  COUPLING_START_TIME=$(tr -d '\r\n' < "$STEADY_FILE")
  if [[ -n "$COUPLING_START_TIME" ]]; then
    export BDM_COUPLING_START_TIME="$COUPLING_START_TIME"
    echo "Coupling start time detected: t=${COUPLING_START_TIME}s"
  fi
fi

deltaT=$(read_control_dict_entry "deltaT")
writeInterval=$(read_control_dict_entry "writeInterval")
if [[ -n "$deltaT" && -n "$writeInterval" ]]; then
  if interval_seconds=$(python3 - <<PY 2>/dev/null
delta_t = float("$deltaT")
write_interval = float("$writeInterval")
print(delta_t * write_interval)
PY
  ); then
    unset BDM_AGENT_EXPORT_INTERVAL
    echo "Agent trajectory export disabled"
  fi
fi

unset BDM_AGENT_OUTPUT_DIR
echo "Agent trajectories are disabled"

PRECICE_RUN_DIR="$MODULAR_COUPLING_DIR/precice-run"
mkdir -p "$PRECICE_RUN_DIR"
echo "Using preCICE exchange directory: $PRECICE_RUN_DIR"
# Do NOT wipe precice-run here: OpenFOAM (acceptor) must start first and may
# already have published connection info. Wipe happens in the OF run.sh.

echo "Capturing BioDynaMo output to $LOG_FILE"

# Set simulation mode if specified (remaining args after participant)
if [[ -n "${1:-}" ]]; then
  export BDM_SIMULATION_MODE="$1"
  echo "Using simulation mode: $1"
  shift
fi

# Build before running
bdm build

# Run BioDynaMo command while teeing output to log file
# Pass the run-specific params CSV as the first argument to the simulation.
#
# Note on ParaView state (.pvsm) generation:
# At the end of the simulation, BioDynaMo invokes its bundled `pvbatch` to run
# generate_pv_state.py, which instantiates a vtkRenderView. The bundled
# ParaView is built against vtkXOpenGLRenderWindow (X11), so:
#   - QT_QPA_PLATFORM=offscreen does NOT help (Qt is not used for rendering).
#   - Unsetting DISPLAY makes pvbatch crash ("bad X server connection. DISPLAY=").
#   - Forwarding the user's DISPLAY/XAUTHORITY (XLaunch / ssh -X) often fails
#     with "Authorization required, but no authorization protocol specified"
#     or "Invalid MIT-MAGIC-COOKIE-1 key" inside the subprocess.
# A successful xdpyinfo probe is not sufficient: an X server may be reachable
# while its GLX vendor/context is broken.  Verify that an OpenGL context can be
# created before entrusting the final PVSM generation to a real DISPLAY.  The
# robust headless fallback is a Mesa-backed virtual X server via xvfb-run.
echo "Launching BioDynaMo run"

display_has_usable_glx() {
  local candidate="$1"

  if ! command -v xdpyinfo >/dev/null 2>&1 ||
     ! DISPLAY="$candidate" xdpyinfo >/dev/null 2>&1; then
    return 1
  fi

  # glxinfo actually creates a GLX context.  This catches displays that
  # xdpyinfo accepts but ParaView later aborts on with X_GLXCreateContext.
  if command -v glxinfo >/dev/null 2>&1; then
    DISPLAY="$candidate" glxinfo -B >/dev/null 2>&1
    return $?
  fi

  # Best available probe on minimal systems without mesa-utils/glxinfo.
  DISPLAY="$candidate" xdpyinfo -queryExtensions 2>/dev/null |
    grep -qE '^[[:space:]]+GLX '
}

find_working_display() {
  local candidate

  if [[ -n "${DISPLAY:-}" ]]; then
    if display_has_usable_glx "$DISPLAY"; then
      printf '%s\n' "$DISPLAY"
      return 0
    fi
  fi

  if command -v xdpyinfo >/dev/null 2>&1; then
    for candidate in :1 :0 :2 :99; do
      if display_has_usable_glx "$candidate"; then
        printf '%s\n' "$candidate"
        return 0
      fi
    done
  fi

  return 1
}

rerun_pvsm_generation() {
  local pvbatch_bin=""
  local pv_state_script=""
  local sim_info_json=""
  local chosen_display=""

  sim_info_json="$CASE_OUTPUT_DIR/cells/simulation_info.json"
  pv_state_script="$BDMSYS/include/core/visualization/paraview/generate_pv_state.py"

  if [[ ! -f "$sim_info_json" || ! -f "$pv_state_script" ]]; then
    return 1
  fi

  if [[ -x "$BDMSYS/third_party/paraview/bin/pvbatch" ]]; then
    pvbatch_bin="$BDMSYS/third_party/paraview/bin/pvbatch"
  elif command -v pvbatch >/dev/null 2>&1; then
    pvbatch_bin="$(command -v pvbatch)"
  fi

  if [[ -z "$pvbatch_bin" ]]; then
    echo "PVSM fallback: no pvbatch executable found"
    return 1
  fi

  if chosen_display="$(find_working_display)"; then
    echo "PVSM fallback: retrying with DISPLAY=$chosen_display via $pvbatch_bin"
    DISPLAY="$chosen_display" "$pvbatch_bin" "$pv_state_script" "$sim_info_json" >> "$LOG_FILE" 2>&1
    return $?
  fi

  if command -v xvfb-run >/dev/null 2>&1; then
    echo "PVSM fallback: no working DISPLAY; retrying with xvfb-run via $pvbatch_bin"
    xvfb-run -a --server-args="-screen 0 1280x1024x24 +extension GLX +render -noreset" \
      "$pvbatch_bin" "$pv_state_script" "$sim_info_json" >> "$LOG_FILE" 2>&1
    return $?
  fi

  echo "PVSM fallback: no working DISPLAY and xvfb-run unavailable"
  return 1
}

BDM_LAUNCH_PREFIX=()
if [[ -n "${BDM_CPUSET:-}" ]]; then
  if command -v taskset >/dev/null 2>&1; then
    BDM_LAUNCH_PREFIX=(taskset -c "$BDM_CPUSET")
    echo "Pinning BioDynaMo participant to cores: $BDM_CPUSET"
  else
    echo "taskset not available; running BioDynaMo without explicit CPU pinning"
  fi
fi

if DISPLAY_WORKING="$(find_working_display)"; then
  echo "Using detected working DISPLAY=$DISPLAY_WORKING for BioDynaMo run"
  DISPLAY="$DISPLAY_WORKING" "${BDM_LAUNCH_PREFIX[@]}" bdm run "$RUN_PARAMS_CONFIG" "$@" 2>&1 | tee "$LOG_FILE"
  run_status="${PIPESTATUS[0]}"
elif command -v xvfb-run >/dev/null 2>&1; then
  echo "No DISPLAY with working GLX detected; using xvfb-run for headless ParaView state generation"
  "${BDM_LAUNCH_PREFIX[@]}" xvfb-run -a \
    --server-args="-screen 0 1280x1024x24 +extension GLX +render -noreset" \
    bdm run "$RUN_PARAMS_CONFIG" "$@" 2>&1 | tee "$LOG_FILE"
  run_status="${PIPESTATUS[0]}"
else
  echo "xvfb-run not found; using current DISPLAY=${DISPLAY:-<unset>} for ParaView state generation"
  echo "  (install xvfb to avoid relying on the user's X server: sudo apt-get install xvfb)"
  "${BDM_LAUNCH_PREFIX[@]}" bdm run "$RUN_PARAMS_CONFIG" "$@" 2>&1 | tee "$LOG_FILE"
  run_status="${PIPESTATUS[0]}"
fi

coupling_completed_before_pvsm() {
  grep -q "Reached end at: final time-window:" "$LOG_FILE" &&
    grep -q "Execution finished after .* simulation_time=" "$LOG_FILE" &&
    grep -q "Close communication channels" "$LOG_FILE"
}

if [[ "$run_status" -ne 0 ]] &&
   grep -q "Error in <ParaviewAdaptor::GenerateParaviewState>" "$LOG_FILE" &&
   coupling_completed_before_pvsm; then
  echo "The coupled simulation completed; only PVSM post-processing failed. Attempting recovery..."
  if rerun_pvsm_generation; then
    echo "PVSM recovery succeeded; overriding run status to success."
    run_status=0
  else
    echo "WARNING: PVSM recovery failed. VTK results are intact; treating the completed coupling as successful."
    run_status=0
  fi
fi

collect_precice_artifacts

echo ""
echo "Run artifacts:"
echo "  Visualization: $CASE_OUTPUT_DIR"
echo "  Trajectories:  disabled"
echo "  Log file:      $LOG_FILE"
echo "  preCICE logs:  $CASE_LOG_DIR"
echo "  preCICE prof.: $CASE_PRECICE_PROFILE_DIR"
echo ""
echo "View commands:"
echo "  ./view_case.sh --list"
echo "  ./view_case.sh $PARTICIPANT"

exit "$run_status"
# bdm view
