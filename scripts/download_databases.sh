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
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/config.sh"
load_config "$CONFIG"

log() { echo "[$(date '+%F %T')] $*"; }

activate_tools

# ----------------------------- Swiss-Prot ------------------------------------
mkdir -p "$DB_DIR"
if [[ ! -s "$DB_DIR/uniprot_sprot.dmnd" || ! -s "$DB_DIR/uniprot_sprot.fasta" ]]; then
    log "Downloading UniProtKB/Swiss-Prot"
    wget -q -O "$DB_DIR/uniprot_sprot.fasta.gz" \
        https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
    gunzip -f "$DB_DIR/uniprot_sprot.fasta.gz"
    wget -q -O "$DB_DIR/reldate.txt" \
        https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/reldate.txt || true
    : # environment activated above
    diamond makedb --in "$DB_DIR/uniprot_sprot.fasta" -d "$DB_DIR/uniprot_sprot"
:
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

# ----------------------------- Reference proteome (optional) -----------------
if [[ -n "${REF_SPECIES:-}" && -z "${REF_PROTEOME:-}" ]]; then
    mkdir -p "$REF_DIR"
    REF_FILE="$REF_DIR/${REF_SPECIES}.pep.all.fa.gz"
    if [[ ! -s "$REF_FILE" ]]; then
        log "Downloading the Ensembl proteome of $REF_SPECIES"
        TMPD=$(mktemp -d "$REF_DIR/dl.XXXX")
        wget -q -r -l1 -nd -np -e robots=off -A '*.pep.all.fa.gz' -P "$TMPD" \
            "https://ftp.ensembl.org/pub/current_fasta/${REF_SPECIES}/pep/"
        SRC=$(ls "$TMPD"/*.pep.all.fa.gz 2>/dev/null | head -1)
        [[ -n "$SRC" ]] || { echo "No proteome found for '$REF_SPECIES' on Ensembl. Check the name at https://ftp.ensembl.org/pub/current_fasta/" >&2; exit 1; }
        echo "$(basename "$SRC")  downloaded $(date '+%F') from Ensembl current_fasta" > "$REF_DIR/${REF_SPECIES}.source.txt"
        mv "$SRC" "$REF_FILE"; rm -rf "$TMPD"
    else
        log "Reference proteome already present, skipping"
    fi
    ls -lh "$REF_FILE"
fi

log "Databases ready:"
ls -lh "$DB_DIR"/uniprot_sprot.{fasta,dmnd} \
       "$EGGNOG_DB"/{eggnog.db,eggnog.taxa.db,eggnog_proteins.dmnd}
