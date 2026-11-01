# Snakemake workflow

The root `Snakefile` is a native Snakemake 9 workflow. It does not wrap the whole
Bash pipeline in one job. Swiss-Prot, InterProScan, eggNOG and the optional
reference search have independent rules, dependencies, resource requests and logs.
Snakemake schedules available branches concurrently and tracks code, parameters,
inputs and incomplete jobs. The original Bash entry point remains available.

## 1. Install on Jarvis

After pushing these files from your Mac and pulling them on Jarvis:

```bash
cd /jarvis/scratch/usr/biello/software/functional-annotation
conda env create -f envs/snakemake.yml
conda activate functional-snakemake
```

This installs the workflow engine and SLURM executor. `--sdm conda` in the commands
below creates a separate annotation environment from `envs/funannot.yml` for the
rules. Create environments on a node with network access. InterProScan and its
member databases must be installed separately; put its **absolute executable path**
in the YAML. Do not point to a wrapper that depends on an unrelated active conda
environment. The annotation environment supplies Java 17. The earlier Bash setup
scripts still use `.sh` configs, not this YAML.

Database downloading/indexing is a separate preparation step: set existing
Swiss-Prot FASTA/index and eggNOG database paths. Avoid modifying shared database
files during jobs. If a DIAMOND index is incompatible, rebuild it separately.

## 2. Configure a project/assembly

```bash
mkdir -p /jarvis/scratch/usr/biello/S_serpentarius/annotation/functional
cp config/snakemake.example.yaml /jarvis/scratch/usr/biello/S_serpentarius/annotation/functional/config.yaml
nano /jarvis/scratch/usr/biello/S_serpentarius/annotation/functional/config.yaml
```

Fill in `genome`, `gff`, `swissprot_fasta`, `swissprot_db`, `eggnog_db` and
`interproscan.command`. The example uses a fresh `functional/snakemake` output
subdirectory so the first run cannot confuse old Bash results with Snakemake outputs.

One YAML config describes one project/species/assembly. For another species or
assembly, copy the config, change the input paths and give it a separate `outdir`
and working directory. If `outdir` is omitted, the default is
`results/<project>/<species>` within the working directory. Unlike the Bash loader,
**relative input/output paths are relative to `--directory`**, not the config file.
Absolute paths are simplest for cluster use. There is no organism list hardcoded
in the rules; the same gene/mRNA/CDS annotation assumptions as the Bash pipeline apply.

## 3. Preview the jobs

```bash
snakemake \
  --snakefile /jarvis/scratch/usr/biello/software/functional-annotation/Snakefile \
  --configfile /jarvis/scratch/usr/biello/S_serpentarius/annotation/functional/config.yaml \
  --directory /jarvis/scratch/usr/biello/S_serpentarius/annotation/functional \
  --cores 16 --dry-run
```

The preview validates the dependency graph and input availability. It does not
validate database contents, executable availability, or cluster submission.

## 4. Run through SLURM

The supplied profile is generic: add `slurm_account` and `slurm_partition` if Jarvis
requires them, following the commented examples in `profiles/slurm/config.yaml`.
Memory, runtime (minutes), and threads are configurable in the project YAML.

```bash
snakemake \
  --snakefile /jarvis/scratch/usr/biello/software/functional-annotation/Snakefile \
  --configfile /jarvis/scratch/usr/biello/S_serpentarius/annotation/functional/config.yaml \
  --directory /jarvis/scratch/usr/biello/S_serpentarius/annotation/functional \
  --profile /jarvis/scratch/usr/biello/software/functional-annotation/profiles/slurm \
  --sdm conda
```

Run the controller in a persistent session or a suitable controller allocation,
according to Jarvis policy. It submits one job per rule. Do not submit this command
as sixteen compute CPUs merely to get parallel workers: the profile requests each
rule's own CPUs/memory. The controller must remain running, and the same repository,
working directory, environments and data paths must be accessible to all workers.

For a local run inside an appropriately sized allocation, omit `--profile` and add
`--cores 16 --resources mem_mb=64000 --rerun-incomplete --sdm conda`.
To use an already prepared software environment instead, omit `--sdm conda` and
ensure every required executable is on PATH on both controller and workers.

Repeat the same command to resume. Snakemake handles incomplete outputs and reruns
rules whose tracked inputs/parameters/code change. A postprocessing-only policy
change does not redo the expensive searches. Keep `.snakemake/` in the working
directory: it stores execution metadata. If you replace external software or
InterProScan member data in place, Snakemake cannot detect that automatically;
record the versions and use `--forcerun interproscan` (or the affected rule).

## Outputs and interpretation

Numbered directories contain longest-isoform GFF/proteins, Swiss-Prot hits,
InterProScan TSV, eggNOG annotations, optional reference RBHs, and final GFF/summary.
Each rule writes a log under `outdir/logs`. This version retains InterProScan TSV;
it does not generate the redundant XML/GFF3 reports. eggNOG intermediate search
files and reference indexes are temporary; final annotation and RBH tables remain.

The final directory also includes the effective YAML configuration and input file
metadata in `*.run.json`. This is not a complete database/software version audit;
record database releases and InterProScan versions alongside the config. A reference
name comparison table is always produced (header only when reference analysis is off).

The final GFF contains the selected longest isoforms and the same annotation types
as the Bash version. eggNOG GO/KEGG/EC/COG annotations remain in the separate eggNOG
table; only the configured fallback names are transferred. Reference RBHs indicate
candidate correspondences, not a phylogenetic demonstration of orthology.

## Tests

In an environment containing Snakemake 9:

```bash
python tests/test_snakemake.py
python tests/test_pipeline.py
```

The Snakemake integration test runs the real workflow engine with fake annotation
executables. It checks the DAG, interrupted-rule recovery, empty-hit outputs,
no-op resume, selective reruns after parameter changes, optional reference analysis,
project isolation, and paths containing spaces. It does not replace a small pilot
run with real annotations/databases on Jarvis. SLURM submission and conda solving
must be validated on that Linux cluster.

References: [Snakemake configuration](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html),
[command-line options](https://snakemake.readthedocs.io/en/stable/executing/cli.html),
and [SLURM executor](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html).
