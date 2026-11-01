# =============================================================================
# Functional annotation pipeline: configuration file
# Copy this file, edit it for your genome, and pass it to the pipeline:
#   cp config/config.example.sh config/my_species.sh
#   sbatch functional_annotation.sh config/my_species.sh
# =============================================================================

# ----------------------------- Inputs ----------------------------------------
PROJECT="my_project"
OUTPUT_ROOT="../results"                  # relative to this config; projects stay separate
SPECIES="mySpecies1"                       # prefix for all output files
GENOME="/path/to/genome.fa"                # assembly used for gene prediction
GFF="/path/to/annotation.gff3"             # structural annotation (BRAKER, AUGUSTUS, ...)
# OUTDIR="/path/to/output"                   # created if missing

# ----------------------------- Databases -------------------------------------
# Folders are created and filled by scripts/download_databases.sh
DB_DIR="/path/to/databases/swissprot"      # uniprot_sprot.fasta + .dmnd
EGGNOG_DB="/path/to/databases/eggnog"      # eggnog.db, eggnog.taxa.db, eggnog_proteins.dmnd

# ----------------------------- Software --------------------------------------
# Command run before conda is used (leave empty if conda is already on PATH)
CONDA_SETUP=""            # e.g. "" or "module load miniconda3"
CONDA_ENV="funannot"                       # envs/funannot.yml: AGAT, DIAMOND, eggNOG-mapper, Java

# InterProScan, installed into the same env by scripts/install_interproscan.sh
IPR_VERSION="5.75-106.0"
IPR_INSTALL_DIR="/path/to/software"        # ~20 GB; unpacked as interproscan-<version>/
IPR_CMD="interproscan.sh"
# Optional: use a cluster module instead of the conda install, e.g.
#   IPR_SETUP="module load openjdk/17.0.8 interproscan/5.75-106.0"
IPR_SETUP=""

# ----------------------------- Resources -------------------------------------
THREADS=${SLURM_CPUS_PER_TASK:-16}

# ----------------------------- InterProScan ----------------------------------
# "yes" if compute nodes have NO internet: disables the EBI precalculated
# match lookup (much slower, but works offline)
IPR_OFFLINE="no"
# Member databases to run; empty = all databases available in your install
IPR_APPL=""                                # e.g. "Pfam,PANTHER,Gene3D,SUPERFAMILY,CDD,SMART"

# ----------------------------- DIAMOND vs Swiss-Prot -------------------------
EVALUE="1e-5"
QCOV=50                                    # min % query coverage
SCOV=50                                    # min % subject coverage
DIAMOND_MODE="--sensitive"                 # --fast | --sensitive | --more-sensitive

# ----------------------------- eggNOG-mapper ---------------------------------
EMAPPER_TAX_SCOPE="auto"                   # e.g. auto | Vertebrata | Aves | Mammalia
EMAPPER_SENSMODE="more-sensitive"          # fast | sensitive | more-sensitive

# ----------------------------- Reference species (optional) ------------------
# Reciprocal best hits against a closely related, well-annotated vertebrate
# from Ensembl. Gives a second, lineage-matched name for each gene and a
# table of where it agrees with the Swiss-Prot name. Leave empty to skip.
# Use the Ensembl directory name (https://ftp.ensembl.org/pub/current_fasta/):
#   birds        gallus_gallus, taeniopygia_guttata
#   mammals      homo_sapiens, mus_musculus
#   reptiles     anolis_carolinensis, chrysemys_picta_bellii
#   amphibians   xenopus_tropicalis
#   fish         danio_rerio, oryzias_latipes
REF_SPECIES=""                             # e.g. "gallus_gallus"
REF_DIR="/path/to/databases/reference"     # proteome saved as <REF_SPECIES>.pep.all.fa.gz
# REF_PROTEOME="/path/to/own.pep.all.fa.gz" # optional: your own Ensembl-format file
# How reference symbols are used in the final GFF3:
#   annotate = add ref_ortholog/ref_symbol/ref_product, keep the Swiss-Prot Name
#   fill     = as annotate, and name genes that have no Swiss-Prot name
#   prefer   = as annotate, and use the reference symbol as Name whenever there
#              is one (the Swiss-Prot name is kept as sprot_name)
REF_NAME_MODE="annotate"

# ----------------------------- Post-processing -------------------------------
# Gene-symbol case. Swiss-Prot best hits are often mouse/rat entries, whose
# symbols are capitalised (Lman1); most vertebrate nomenclatures use uppercase.
#   upper = LMAN1 (birds, reptiles, amphibians, human/non-rodent mammals)
#   keep  = leave symbols as they come from Swiss-Prot
#   (for mouse/rat genomes use "keep"; for zebrafish/fish see README)
NAME_CASE="keep"
# Give a name from eggNOG (Preferred_name) to genes without a Swiss-Prot name
EGGNOG_NAME_FALLBACK="yes"
# Remove the "(InterPro)"/"(PANTHER)" source tags from GO terms, which strict
# GFF3 validators and NCBI submission reject
CLEAN_GO_TAGS="yes"
