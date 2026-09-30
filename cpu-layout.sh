#!/bin/bash

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

build_physical_cpuset_from_slice() {
  local start="${1:-0}"
  local length="${2:-1}"
  local -a cpu_list
  local -a slice
  local IFS=,

  mapfile -t cpu_list < <(detect_physical_cpu_list)
  slice=("${cpu_list[@]:start:length}")
  echo "${slice[*]}"
}

build_physical_cpuset_excluding_cpu() {
  local excluded_cpu="${1:-0}"
  local length="${2:-1}"
  local anchor_cpu
  local -a cpu_list
  local -a selected=()
  local cpu
  local IFS=,

  anchor_cpu="$(detect_physical_anchor_cpu_for_cpu "$excluded_cpu")"
  mapfile -t cpu_list < <(detect_physical_cpu_list)

  for cpu in "${cpu_list[@]}"; do
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