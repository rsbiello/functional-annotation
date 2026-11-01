#!/bin/bash
# Run several configs locally, or submit one independent SLURM job per config.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MODE=${1:?Usage: run_projects.sh local|slurm config.sh [config.sh ...]}
shift
[[ $# -gt 0 && "$MODE" =~ ^(local|slurm)$ ]] || { echo 'Provide local|slurm and at least one config' >&2; exit 1; }
for input in "$@"; do
    config="$(cd "$(dirname "$input")" && pwd)/$(basename "$input")"
    if [[ "$MODE" == local ]]; then
        bash "$HERE/functional_annotation.sh" "$config"
    else
        # --chdir lets the copied SLURM script locate its helper scripts.
        sbatch --chdir="$HERE" "$HERE/functional_annotation.sh" "$config"
    fi
done
