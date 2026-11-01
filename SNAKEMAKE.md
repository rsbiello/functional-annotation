# Run functional annotation with Snakemake on Jarvis

Use **one Conda environment: `funannot`**. It contains Snakemake, the SLURM executor,
AGAT, DIAMOND, eggNOG-mapper, Java and Python. Activate it for both setup and runs.
InterProScan and the databases are installed in shared folders outside Conda and
reused for different species. The workflow uses the active environment; it does
not create per-rule environments. Do not add `--sdm conda` to these commands.

## 1. Install or update the environment

After pushing the repository changes from your Mac, run on Jarvis:

```bash
cd /jarvis/scratch/usr/biello/software/functional-annotation
git pull --ff-only
conda env create -f envs/funannot.yml
conda activate funannot
```

If `funannot` already exists, replace the creation command with:

```bash
conda env update -n funannot -f envs/funannot.yml
conda activate funannot
```

If you previously created `functional-snakemake`, leave it unused; it is not needed
for this setup. No environment or downloaded database needs to be deleted.
Create/update the environment on a suitable node with internet access. All workers
must be able to access this same Conda environment, repository and data directories.

## 2. Install shared InterProScan and databases once

Skip installations you have already completed at these paths. Create the setup
configuration from the repository directory:

```bash
cat > config/setup.sh <<'CONFIG'
PROJECT="shared_databases"
SPECIES="shared"
CONDA_SETUP=""
CONDA_ENV="funannot"
DB_DIR="/jarvis/scratch/usr/biello/databases/swissprot"
EGGNOG_DB="/jarvis/scratch/usr/biello/databases/eggnog"
IPR_VERSION="5.75-106.0"
IPR_INSTALL_DIR="/jarvis/scratch/usr/biello/software"
REF_SPECIES=""
CONFIG

bash scripts/install_interproscan.sh config/setup.sh
bash scripts/download_databases.sh config/setup.sh
```

These scripts use Bash configuration (`setup.sh`). The annotation workflow uses
YAML (`config.yaml`, below). No genome/GFF is required for database setup.
Downloads, extraction and indexing need substantial disk space and a suitable node
with internet access. Stop and resolve any failed command before proceeding.
The installer checks InterProScan's archive checksum and runs a small Pfam test.
Keep eggNOG-mapper 2.1.15 paired with the downloaded 5.0.2 databases.

Verify the tools and Swiss-Prot index:

```bash
snakemake --version
interproscan.sh -version
diamond dbinfo -d /jarvis/scratch/usr/biello/databases/swissprot/uniprot_sprot.dmnd
ls -lh /jarvis/scratch/usr/biello/databases/eggnog/
```

The eggNOG directory must contain `eggnog.db`, `eggnog.taxa.db`, and
`eggnog_proteins.dmnd`. Do not modify shared databases while annotation jobs run.

## 3. Configure the species/assembly

```bash
mkdir -p /jarvis/scratch/usr/biello/S_serpentarius/annotation/functional
cp config/snakemake.example.yaml /jarvis/scratch/usr/biello/S_serpentarius/annotation/functional/config.yaml
nano /jarvis/scratch/usr/biello/S_serpentarius/annotation/functional/config.yaml
```

If the project config already exists, edit it instead of copying over it. Set its
real genome and GFF paths. Update the existing database/software entries to:

```yaml
swissprot_fasta: /jarvis/scratch/usr/biello/databases/swissprot/uniprot_sprot.fasta
swissprot_db: /jarvis/scratch/usr/biello/databases/swissprot/uniprot_sprot.dmnd
eggnog_db: /jarvis/scratch/usr/biello/databases/eggnog
interproscan:
  command: /jarvis/scratch/usr/biello/software/interproscan-5.75-106.0/interproscan.sh
  offline: true
  applications: ""
```

Keep the other example settings. `offline: true` disables the online match lookup,
allowing annotation on compute nodes without internet access.

Use one YAML per assembly/run. Change `project`, `species`, the inputs and `outdir`
for other projects. Use separate working and output directories for independent
runs. The example writes to `functional/snakemake` to avoid reusing old Bash outputs.
Relative data paths are resolved against `--directory`, not the YAML's location;
absolute paths are simplest. Eukaryotic gene/mRNA/CDS annotations are assumed.

## 4. Preview and run

Activate the same environment and define these shortcuts in each new shell:

```bash
conda activate funannot
REPO=/jarvis/scratch/usr/biello/software/functional-annotation
RUN=/jarvis/scratch/usr/biello/S_serpentarius/annotation/functional
```

Preview the jobs:

```bash
snakemake \
  --snakefile "$REPO/Snakefile" \
  --configfile "$RUN/config.yaml" \
  --directory "$RUN" \
  --cores 16 --dry-run
```

If Jarvis requires a SLURM account or partition, edit the commented examples in
`profiles/slurm/config.yaml`. Memory, CPU count and runtime (minutes) are configurable
in the project YAML. Then submit:

```bash
snakemake \
  --snakefile "$REPO/Snakefile" \
  --configfile "$RUN/config.yaml" \
  --directory "$RUN" \
  --profile "$REPO/profiles/slurm"
```

Keep the controller running in a persistent session or controller allocation
permitted by Jarvis policy. It submits one job per rule; searches can run concurrently.
Workers must inherit the activated environment's PATH. The standard SLURM executor
exports the submission environment; site-specific settings must not strip it.

For local execution inside a suitable compute allocation:

```bash
snakemake \
  --snakefile "$REPO/Snakefile" \
  --configfile "$RUN/config.yaml" \
  --directory "$RUN" \
  --cores 16 --resources mem_mb=64000 --rerun-incomplete
```

Repeat the same execution command to resume. Snakemake tracks input, parameter and
code changes, and the profile enables incomplete-job recovery. Keep `.snakemake/`
in the working directory. A postprocessing-only change does not repeat searches.
Because software is managed externally, changing the active environment or replacing
InterProScan/database data in place may require `--forcerun RULE` for the affected
rule, or a fresh output directory. Record software and database versions for each run.

## Outputs and validation

Numbered directories contain longest-isoform GFF/proteins, Swiss-Prot hits,
InterProScan TSV, eggNOG annotations, optional reference RBHs, and final GFF/summary.
Logs are in `outdir/logs`. Final results are in `outdir/06_final`. InterProScan XML/GFF3
and intermediate eggNOG search files are not retained by this workflow.

`*.run.json` records the effective config and input metadata; it is not a complete
software/database version audit. The reference comparison table is header-only
when reference analysis is disabled. eggNOG GO/KEGG/EC/COG annotations remain in its
separate table; only configured fallback names are transferred into the final GFF.
Reference RBHs identify candidate correspondences, not proven orthology.

```bash
python tests/test_snakemake.py
python tests/test_pipeline.py
```

Tests run real Snakemake with simulated annotation executables. They cover dependency
resolution, interrupted jobs, empty/positive hits, resume, selective parameter reruns,
optional references, project isolation and paths with spaces. A small real-data pilot,
Conda solving and SLURM submission still need validation on Jarvis.

References: [Snakemake CLI](https://snakemake.readthedocs.io/en/stable/executing/cli.html)
and [SLURM executor](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html).
