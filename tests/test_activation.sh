#!/bin/bash
set -eo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
source "$ROOT/scripts/config.sh"
CONDA_SETUP=""
CONDA_ENV=funannot
# A preactivated environment must not invoke Conda at all.
(
    set -u
    CONDA_PREFIX=/envs/funannot CONDA_DEFAULT_ENV=funannot
    conda() { return 99; }
    activate_tools
    [[ $- == *u* ]]
)
# Simulate a third-party hook reading an unset variable.
(
    set -u
    unset CONDA_PREFIX CONDA_DEFAULT_ENV xml_catalog_files_libxml2
    conda() {
        if [[ "$1" == shell.bash ]]; then
            printf '%s\n' 'dummy="$xml_catalog_files_libxml2"'
        else
            dummy="$xml_catalog_files_libxml2"
            export CONDA_PREFIX=/envs/funannot
        fi
    }
    activate_tools
    [[ "$CONDA_PREFIX" == /envs/funannot && $- == *u* ]]
)
# Failure must propagate, and strict mode must be restored.
(
    set -u
    unset CONDA_PREFIX CONDA_DEFAULT_ENV
    conda() { if [[ "$1" == shell.bash ]]; then echo ':'; else return 42; fi; }
    status=0
    activate_tools || status=$?
    [[ "$status" == 42 && $- == *u* ]]
)
# Do not enable nounset if the caller did not request it.
(
    set +u
    unset CONDA_PREFIX CONDA_DEFAULT_ENV
    conda() { if [[ "$1" == shell.bash ]]; then echo ':'; fi; }
    activate_tools
    [[ $- != *u* ]]
)
echo 'PASS: existing environment, unset hook variables, failure propagation, shell option restoration'
