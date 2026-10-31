#!/bin/bash
# =============================================================================
# Install InterProScan into the pipeline's conda environment
#
# Usage:  bash scripts/install_interproscan.sh config/my_species.sh
#
# Downloads the official release from EBI, checks the MD5, unpacks it into
# IPR_INSTALL_DIR, indexes the member databases (setup.py) and puts an
# `interproscan.sh` wrapper in the conda environment's bin/, so that
#   conda activate funannot && interproscan.sh -version
# works with no modules. Java 17 comes from the same environment.
#
# Needs internet (run it on the login node) and ~20 GB of disk.
# InterProScan runs on Linux only.
# =============================================================================
set -euo pipefail

CONFIG=${1:?Usage: bash scripts/install_interproscan.sh <config.sh>}
# shellcheck source=/dev/null
source "$CONFIG"

IPR_VERSION=${IPR_VERSION:-5.75-106.0}
IPR_INSTALL_DIR=${IPR_INSTALL_DIR:?set IPR_INSTALL_DIR in the config}

log() { echo "[$(date '+%F %T')] $*"; }

[[ -n "${CONDA_SETUP:-}" ]] && eval "$CONDA_SETUP"
eval "$(conda shell.bash hook)"
conda activate "$CONDA_ENV"

java -version 2>&1 | head -1

TARBALL="interproscan-${IPR_VERSION}-64-bit.tar.gz"
URL="https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/${IPR_VERSION}"
IPR_HOME="$IPR_INSTALL_DIR/interproscan-${IPR_VERSION}"

mkdir -p "$IPR_INSTALL_DIR"
cd "$IPR_INSTALL_DIR"

if [[ ! -x "$IPR_HOME/interproscan.sh" ]]; then
    log "Downloading InterProScan $IPR_VERSION (several GB)"
    wget -c "$URL/$TARBALL"
    wget -c "$URL/$TARBALL.md5"
    md5sum -c "$TARBALL.md5"

    log "Unpacking"
    tar -pxzf "$TARBALL"
    rm -f "$TARBALL"

    log "Indexing member databases (setup.py)"
    cd "$IPR_HOME"
    python3 setup.py -f interproscan.properties
else
    log "InterProScan already unpacked in $IPR_HOME"
fi

# Wrapper inside the conda environment
cat > "$CONDA_PREFIX/bin/interproscan.sh" <<EOF
#!/bin/bash
exec "$IPR_HOME/interproscan.sh" "\$@"
EOF
chmod +x "$CONDA_PREFIX/bin/interproscan.sh"
log "Wrapper written to $CONDA_PREFIX/bin/interproscan.sh"

log "Quick test on the bundled test proteins"
cd "$IPR_HOME"
interproscan.sh -i test_all_appl.fasta -f tsv -appl Pfam -dp -o "$IPR_INSTALL_DIR/ipr_test.tsv" >/dev/null
log "Test OK: $(wc -l < "$IPR_INSTALL_DIR/ipr_test.tsv") Pfam matches"
interproscan.sh -version | head -3
