"""Functional annotation: one isolated project/assembly per configuration."""
from pathlib import Path
import re
from snakemake.utils import min_version

min_version("9.0")
ROOT = Path(workflow.basedir).resolve()
# Relative data paths follow Snakemake's working directory (--directory).
for key in ("project", "species", "genome", "gff", "swissprot_fasta", "swissprot_db", "eggnog_db"):
    if not config.get(key):
        raise ValueError(f"Missing required configuration key: {key}")
for key in ("project", "species"):
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", str(config[key])):
        raise ValueError(f"Invalid {key} identifier")
P = config["species"]
OUT = str(Path(config.get("outdir", f"results/{config['project']}/{P}")))
def out(relative):
    return str(Path(OUT) / relative)
def helper(name):
    return str(ROOT / "scripts" / name)
ENV = str(ROOT / "envs" / "funannot.yml")
SEARCH = config.get("search", {})
IPR = config.get("interproscan", {})
EGG = config.get("eggnog", {})
REF = config.get("reference", {}) or {}
POST = config.get("postprocess", {})
if POST.get("name_case", "keep") not in ("keep", "upper"):
    raise ValueError("postprocess.name_case must be keep or upper")
if REF.get("name_mode", "annotate") not in ("annotate", "fill", "prefer"):
    raise ValueError("reference.name_mode must be annotate, fill or prefer")
if REF and (not REF.get("proteome") or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", REF.get("species", ""))):
    raise ValueError("reference needs a proteome path and a safe species identifier")
for section, key in ((IPR, "offline"), (POST, "eggnog_names"), (POST, "clean_go_tags")):
    if key in section and not isinstance(section[key], bool):
        raise ValueError(f"{key} must be a YAML boolean (true or false)")
MODE = SEARCH.get("mode", "sensitive")
if MODE not in ("fast", "mid-sensitive", "sensitive", "more-sensitive", "very-sensitive", "ultra-sensitive"):
    raise ValueError("Invalid search.mode")
EVALUE = float(SEARCH.get("evalue", 1e-5))
QCOV = float(SEARCH.get("query_cover", 50))
SCOV = float(SEARCH.get("subject_cover", 50))
if EVALUE <= 0 or not (0 <= QCOV <= 100 and 0 <= SCOV <= 100):
    raise ValueError("Require positive evalue and coverage percentages between 0 and 100")
NTHREADS = int(config.get("threads", 16))
if NTHREADS < 1:
    raise ValueError("threads must be positive")
MEMORY = config.get("memory_mb", {})
RUNTIME = config.get("runtime_minutes", {})
def mem(step, default):
    return int(MEMORY.get(step, default))
def runtime(step, default):
    return int(RUNTIME.get(step, default))

LONGEST = out(f"01_proteins/{P}.longest.gff3")
PROTEINS = out(f"01_proteins/{P}.longest.prot.fa")
SPROT = out(f"02_swissprot/{P}.sprot.outfmt6.tsv")
IPRTSV = out(f"03_interproscan/{P}.ipr.tsv")
EGGANN = out(f"04_eggnog/{P}.emapper.annotations")
MERGED = out(f"05_merged/{P}.gff3")
RBH = out("05_reference/rbh.tsv")
FINAL = out(f"06_final/{P}.functional.gff3")
SUMMARY = out(f"06_final/{P}.functional_annotation_summary.tsv")
COMPARISON = out(f"06_final/{P}.reference_name_check.tsv")
PROVENANCE = out(f"06_final/{P}.run.json")

rule all:
    input:
        FINAL, SUMMARY, COMPARISON, PROVENANCE

rule longest_isoform:
    input:
        config["gff"]
    output:
        LONGEST
    log:
        out("logs/longest_isoform.log")
    resources:
        mem_mb=mem("agat", 16000), runtime=runtime("agat", 120)
    conda:
        ENV
    shell:
        "agat_sp_keep_longest_isoform.pl --gff {input:q} --output {output:q} > {log:q} 2>&1"

rule extract_proteins:
    input:
        gff=LONGEST, genome=config["genome"], cleaner=helper("clean_proteins.py")
    output:
        PROTEINS
    log:
        out("logs/extract_proteins.log")
    resources:
        mem_mb=mem("agat", 16000), runtime=runtime("agat", 120)
    conda:
        ENV
    shell:
        r"""
        scratch=$(mktemp -d)
        trap 'rm -rf "$scratch"' EXIT
        agat_sp_extract_sequences.pl --gff {input.gff:q} --fasta {input.genome:q} --protein --output "$scratch/raw.fa" > {log:q} 2>&1
        python3 {input.cleaner:q} "$scratch/raw.fa" {output:q} >> {log:q} 2>&1
        """

rule swissprot:
    input:
        proteins=PROTEINS, db=config["swissprot_db"]
    output:
        SPROT
    params:
        evalue=EVALUE, qcov=QCOV, scov=SCOV, mode="--" + MODE
    threads:
        NTHREADS
    resources:
        mem_mb=mem("swissprot", 16000), runtime=runtime("swissprot", 240)
    log:
        out("logs/swissprot.log")
    conda:
        ENV
    shell:
        r"""
        scratch=$(mktemp -d)
        trap 'rm -rf "$scratch"' EXIT
        diamond dbinfo -d {input.db:q} > {log:q} 2>&1
        diamond blastp -q {input.proteins:q} -d {input.db:q} -o {output:q} --outfmt 6 \
          --evalue {params.evalue} --query-cover {params.qcov} --subject-cover {params.scov} \
          --max-target-seqs 1 {params.mode:q} -p {threads} --tmpdir "$scratch" >> {log:q} 2>&1
        """

rule interproscan:
    input:
        PROTEINS
    output:
        IPRTSV
    params:
        command=IPR.get("command", "interproscan.sh"),
        offline="-dp" if IPR.get("offline", False) else "",
        applications=(["-appl", IPR["applications"]] if IPR.get("applications") else [])
    threads:
        NTHREADS
    resources:
        mem_mb=mem("interproscan", 64000), runtime=runtime("interproscan", 4320)
    log:
        out("logs/interproscan.log")
    conda:
        ENV
    shell:
        r"""
        scratch=$(mktemp -d)
        trap 'rm -rf "$scratch"' EXIT
        {params.command:q} -i {input:q} -o {output:q} -f TSV -iprlookup -goterms \
          -cpu {threads} {params.offline} {params.applications:q} -T "$scratch" > {log:q} 2>&1
        """

rule eggnog:
    input:
        proteins=PROTEINS,
        db=str(Path(config["eggnog_db"]) / "eggnog.db"),
        taxa=str(Path(config["eggnog_db"]) / "eggnog.taxa.db"),
        diamond=str(Path(config["eggnog_db"]) / "eggnog_proteins.dmnd")
    output:
        EGGANN
    params:
        database=config["eggnog_db"], tax_scope=str(EGG.get("tax_scope", "auto")),
        sensitivity=EGG.get("sensitivity", "more-sensitive"), prefix=P
    threads:
        NTHREADS
    resources:
        mem_mb=mem("eggnog", 64000), runtime=runtime("eggnog", 1440)
    log:
        out("logs/eggnog.log")
    conda:
        ENV
    shell:
        r"""
        scratch=$(mktemp -d)
        trap 'rm -rf "$scratch"' EXIT
        emapper.py -i {input.proteins:q} --itype proteins -m diamond --sensmode {params.sensitivity:q} \
          --data_dir {params.database:q} --tax_scope {params.tax_scope:q} -o {params.prefix:q} \
          --output_dir "$scratch" --temp_dir "$scratch" --cpu {threads} --override > {log:q} 2>&1
        cp "$scratch/"{params.prefix:q}".emapper.annotations" {output:q}
        """

rule merge:
    input:
        gff=LONGEST, blast=SPROT, ipr=IPRTSV, fasta=config["swissprot_fasta"]
    output:
        MERGED
    params:
        evalue=EVALUE
    resources:
        mem_mb=mem("agat", 16000), runtime=runtime("agat", 120)
    log:
        out("logs/merge.log")
    conda:
        ENV
    shell:
        r"""
        scratch=$(mktemp -d)
        trap 'rm -rf "$scratch"' EXIT
        agat_sp_manage_functional_annotation.pl --gff {input.gff:q} --blast {input.blast:q} \
          --db {input.fasta:q} --blast_evalue {params.evalue} --interpro {input.ipr:q} \
          --output "$scratch/agat" > {log:q} 2>&1
        shopt -s nullglob
        candidates=("$scratch/agat/"*.gff*)
        if [[ ${{#candidates[@]}} -ne 1 ]]; then
            echo 'Expected exactly one AGAT output GFF' >> {log:q}; exit 1
        fi
        cp "${{candidates[0]}}" {output:q}
        """

if REF:
    rule reference_rbh:
        input:
            proteins=PROTEINS, reference=REF["proteome"],
            runner=helper("reference_rbh.sh"), parser=helper("prepare_reference.awk")
        output:
            RBH
        params:
            evalue=EVALUE, qcov=QCOV, scov=SCOV
        threads:
            NTHREADS
        resources:
            mem_mb=mem("reference", 16000), runtime=runtime("reference", 240)
        log:
            out("logs/reference_rbh.log")
        conda:
            ENV
        shell:
            r"""
            scratch=$(mktemp -d)
            trap 'rm -rf "$scratch"' EXIT
            bash {input.runner:q} {input.proteins:q} {input.reference:q} "$scratch" \
              {threads} {params.evalue} {params.qcov} {params.scov} > {log:q} 2>&1
            cp "$scratch/rbh.tsv" {output:q}
            """

rule postprocess:
    input:
        gff=MERGED, eggnog=EGGANN, rbh=([RBH] if REF else []),
        reference_awk=helper("apply_reference_names.awk"), eggnog_awk=helper("add_eggnog_names.awk")
    output:
        gff=FINAL, comparison=COMPARISON
    params:
        name_case=POST.get("name_case", "keep"),
        eggnog_names=POST.get("eggnog_names", True), clean_go_tags=POST.get("clean_go_tags", True),
        reference_species=REF.get("species", ""), reference_mode=REF.get("name_mode", "annotate")
    resources:
        mem_mb=mem("postprocess", 8000), runtime=runtime("postprocess", 60)
    log:
        out("logs/postprocess.log")
    conda:
        ENV
    script:
        "workflow/scripts/postprocess.py"

rule summary:
    input:
        proteins=PROTEINS, blast=SPROT, ipr=IPRTSV, eggnog=EGGANN, gff=FINAL,
        comparison=COMPARISON, rbh=([RBH] if REF else [])
    output:
        summary=SUMMARY, provenance=PROVENANCE
    params:
        settings=config
    resources:
        mem_mb=mem("postprocess", 8000), runtime=runtime("postprocess", 60)
    log:
        out("logs/summary.log")
    conda:
        ENV
    script:
        "workflow/scripts/summary.py"
