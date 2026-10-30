# =============================================================================
# Functional annotation pipeline: configuration file
# Copy this file, edit it for your genome, and pass it to the pipeline:
#   cp config/config.example.sh config/my_species.sh
#   sbatch functional_annotation.sh config/my_species.sh
# =============================================================================

# ----------------------------- Inputs ----------------------------------------
SPECIES="mySpecies1"                       # prefix for all output files
GENOME="/path/to/genome.fa"                # assembly used for gene prediction
GFF="/path/to/annotation.gff3"             # structural annotation (BRAKER, AUGUSTUS, ...)
OUTDIR="/path/to/output"                   # created if missing

# ----------------------------- Databases -------------------------------------
# Folders are created and filled by scripts/download_databases.sh
DB_DIR="/path/to/databases/swissprot"      # uniprot_sprot.fasta + .dmnd
EGGNOG_DB="/path/to/databases/eggnog"      # eggnog.db, eggnog.taxa.db, eggnog_proteins.dmnd

# ----------------------------- Software --------------------------------------
# Commands run before the conda environments are activated (leave empty if
# conda is already on your PATH)
CONDA_SETUP="module load conda"            # e.g. "" or "module load miniconda3"
CONDA_ANNOT="annot"                        # env with diamond + agat (envs/annot.yml)
CONDA_EMAPPER="emapper"                    # env with eggnog-mapper (envs/emapper.yml)

# Commands that make interproscan.sh available (module system, PATH, ...)
IPR_SETUP="module load openjdk interproscan"
# e.g. IPR_SETUP="module load openjdk/17.0.8 interproscan/5.75-106.0"
# e.g. IPR_SETUP="export PATH=/opt/interproscan-5.75-106.0:\$PATH"
IPR_CMD="interproscan.sh"

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

# ----------------------------- Post-processing -------------------------------
# Gene-symbol case. Swiss-Prot best hits are often mouse/rat entries, whose
# symbols are capitalised (Lman1); most vertebrate nomenclatures use uppercase.
#   upper = LMAN1 (birds, reptiles, amphibians, human/non-rodent mammals)
#   keep  = leave symbols as they come from Swiss-Prot
#   (for mouse/rat genomes use "keep"; for zebrafish/fish see README)
NAME_CASE="upper"
# Give a name from eggNOG (Preferred_name) to genes without a Swiss-Prot name
EGGNOG_NAME_FALLBACK="yes"
# Remove the "(InterPro)"/"(PANTHER)" source tags from GO terms, which strict
# GFF3 validators and NCBI submission reject
CLEAN_GO_TAGS="yes"
