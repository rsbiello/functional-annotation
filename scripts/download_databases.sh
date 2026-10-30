#!/bin/bash
# =============================================================================
# Download and index the databases used by the pipeline
#   - UniProtKB/Swiss-Prot (FASTA + DIAMOND index)       < 1 GB
#   - eggNOG-mapper v5.0.2 data                           ~ 50 GB unpacked
#
# Usage:  bash scripts/download_databases.sh config/my_species.sh
#
# Run it on a node with internet access (often the login node). Files that
# already exist are skipped. eggNOG files are fetched directly from
# eggnog5.embl.de, because download_eggnog_data.py in some eggNOG-mapper
# releases still points to the retired eggnogdb.embl.de host.
# =============================================================================
set -euo pipefail

CONFIG=${1:?Usage: bash scripts/download_databases.sh <config.sh>}
# shellcheck source=/dev/null
source "$CONFIG"

log() { echo "[$(date '+%F %T')] $*"; }

[[ -n "${CONDA_SETUP:-}" ]] && eval "$CONDA_SETUP"
eval "$(conda shell.bash hook)"

# ----------------------------- Swiss-Prot ------------------------------------
mkdir -p "$DB_DIR"
if [[ ! -s "$DB_DIR/uniprot_sprot.dmnd" ]]; then
    log "Downloading UniProtKB/Swiss-Prot"
    wget -q -O "$DB_DIR/uniprot_sprot.fasta.gz" \
        https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
    gunzip -f "$DB_DIR/uniprot_sprot.fasta.gz"
    wget -q -O "$DB_DIR/reldate.txt" \
        https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/reldate.txt || true
    conda activate "$CONDA_ANNOT"
    diamond makedb --in "$DB_DIR/uniprot_sprot.fasta" -d "$DB_DIR/uniprot_sprot"
    conda deactivate
else
    log "Swiss-Prot already present, skipping"
fi

# ----------------------------- eggNOG ----------------------------------------
mkdir -p "$EGGNOG_DB"
cd "$EGGNOG_DB"
BASE="http://eggnog5.embl.de/download/emapperdb-5.0.2"

if [[ ! -s eggnog.db ]]; then
    log "Downloading eggnog.db (large, be patient)"
    rm -f eggnog.db.gz
    wget -q "$BASE/eggnog.db.gz" && gunzip eggnog.db.gz
fi
if [[ ! -s eggnog.taxa.db ]]; then
    log "Downloading eggnog.taxa.db"
    rm -f eggnog.taxa.tar.gz
    wget -q "$BASE/eggnog.taxa.tar.gz" && tar -zxf eggnog.taxa.tar.gz && rm eggnog.taxa.tar.gz
fi
if [[ ! -s eggnog_proteins.dmnd ]]; then
    log "Downloading eggnog_proteins.dmnd"
    rm -f eggnog_proteins.dmnd.gz
    wget -q "$BASE/eggnog_proteins.dmnd.gz" && gunzip eggnog_proteins.dmnd.gz
fi

log "Databases ready:"
ls -lh "$DB_DIR"/uniprot_sprot.{fasta,dmnd} \
       "$EGGNOG_DB"/{eggnog.db,eggnog.taxa.db,eggnog_proteins.dmnd}
