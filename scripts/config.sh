# Shared configuration loader. Config files are trusted executable Bash.
load_config() {
    CONFIG="$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
    source "$CONFIG"
    local key value
    for key in GENOME GFF DB_DIR EGGNOG_DB REF_DIR REF_PROTEOME IPR_INSTALL_DIR OUTPUT_ROOT OUTDIR; do
        value=${!key:-}
        if [[ -n "$value" && "$value" != /* ]]; then
            printf -v "$key" '%s/%s' "$(dirname "$CONFIG")" "$value"
        fi
    done
    PROJECT=${PROJECT:-default}
    SPECIES=${SPECIES:?Set SPECIES in config}
    [[ "$PROJECT" =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*$ && "$SPECIES" =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*$ ]] || { echo 'Invalid PROJECT or SPECIES identifier' >&2; return 1; }
    OUTDIR=${OUTDIR:-${OUTPUT_ROOT:-$(dirname "$CONFIG")/results}/$PROJECT/$SPECIES}
    if [[ -n "${REF_SPECIES:-}" ]]; then
        [[ "$REF_SPECIES" =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*$ ]] || { echo 'Invalid REF_SPECIES identifier' >&2; return 1; }
    fi
    if [[ "${IPR_CMD:-}" == */* && "$IPR_CMD" != /* ]]; then
        IPR_CMD="$(dirname "$CONFIG")/$IPR_CMD"
    fi
    THREADS=${THREADS:-${SLURM_CPUS_PER_TASK:-16}}
    [[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || { echo 'THREADS must be a positive integer' >&2; return 1; }
    NAME_CASE=${NAME_CASE:-keep}
    [[ "$NAME_CASE" =~ ^(keep|upper)$ ]] || return 1
    CONDA_ENV=${CONDA_ENV-funannot}
    EVALUE=${EVALUE:-1e-5}; QCOV=${QCOV:-50}; SCOV=${SCOV:-50}
    DIAMOND_MODE=${DIAMOND_MODE:---sensitive}
    [[ "$DIAMOND_MODE" =~ ^--(fast|mid-sensitive|sensitive|more-sensitive|very-sensitive|ultra-sensitive)$ ]] || return 1
    IPR_CMD=${IPR_CMD:-interproscan.sh}; IPR_OFFLINE=${IPR_OFFLINE:-no}
    EMAPPER_TAX_SCOPE=${EMAPPER_TAX_SCOPE:-auto}; EMAPPER_SENSMODE=${EMAPPER_SENSMODE:-more-sensitive}
    EGGNOG_NAME_FALLBACK=${EGGNOG_NAME_FALLBACK:-yes}; CLEAN_GO_TAGS=${CLEAN_GO_TAGS:-yes}
    for key in IPR_OFFLINE EGGNOG_NAME_FALLBACK CLEAN_GO_TAGS; do
        [[ "${!key}" =~ ^(yes|no)$ ]] || { echo "$key must be yes or no" >&2; return 1; }
    done
}
activate_tools() {
    [[ -z "${CONDA_SETUP:-}" ]] || eval "$CONDA_SETUP"
    if [[ -n "$CONDA_ENV" ]]; then
        eval "$(conda shell.bash hook)"
        conda activate "$CONDA_ENV"
    fi
}
