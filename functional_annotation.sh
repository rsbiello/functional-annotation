#!/bin/bash
#SBATCH --job-name=func_annot
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=3-00:00:00
#SBATCH --output=func_annot_%j.out
#SBATCH --error=func_annot_%j.err

# =============================================================================
# Functional annotation of a vertebrate genome
#
#   1. Longest isoform per gene + protein extraction    (AGAT)
#   2. Gene names from UniProtKB/Swiss-Prot              (DIAMOND)
#   3. Protein domains and GO terms                      (InterProScan)
#   4. Orthology-based names, GO, KEGG, EC, COG          (eggNOG-mapper)
#   5. Merge names + InterPro/GO into the GFF3           (AGAT)
#   5b. Optional: orthologs in a reference species       (DIAMOND, reciprocal best hits)
#   6. Post-processing: name case, reference names, eggNOG fallback names,
#      clean GO tags
#   7. Summary statistics and software versions
#
# Usage:
#   sbatch functional_annotation.sh config/my_species.sh     # SLURM
#   bash   functional_annotation.sh config/my_species.sh     # any Linux machine
#
# Each step is skipped when its output exists, so the job can be resubmitted
# after a crash or a timeout and continues where it stopped.
# =============================================================================

set -euo pipefail
export LC_NUMERIC=C          # decimal point, not comma, in the summary

CONFIG=${1:?Usage: functional_annotation.sh <config.sh>}
# shellcheck source=/dev/null
source "$CONFIG"

# Directory of this script (works under SLURM, which copies the script)
if [[ -n "${SLURM_SUBMIT_DIR:-}" && -d "$SLURM_SUBMIT_DIR/scripts" ]]; then
    PIPE_DIR="$SLURM_SUBMIT_DIR"
else
    PIPE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi

log() { echo "[$(date '+%F %T')] $*"; }
die() { echo "ERROR: $*" >&2; exit 1; }

# ----------------------------- Checks ----------------------------------------
[[ -s "$GENOME" ]] || die "genome not found: $GENOME"
[[ -s "$GFF"    ]] || die "GFF3 not found: $GFF"
[[ -s "$DB_DIR/uniprot_sprot.dmnd" ]] \
    || die "Swiss-Prot missing in $DB_DIR. Run: bash scripts/download_databases.sh $CONFIG"
[[ -s "$EGGNOG_DB/eggnog.db" && -s "$EGGNOG_DB/eggnog_proteins.dmnd" ]] \
    || die "eggNOG data missing in $EGGNOG_DB. Run: bash scripts/download_databases.sh $CONFIG"

# Optional reference species for reciprocal best hits (see config)
REF_SPECIES=${REF_SPECIES:-}
REF_NAME_MODE=${REF_NAME_MODE:-annotate}
if [[ -n "$REF_SPECIES" ]]; then
    REF_PROTEOME=${REF_PROTEOME:-$REF_DIR/${REF_SPECIES}.pep.all.fa.gz}
    [[ -s "$REF_PROTEOME" ]] \
        || die "reference proteome missing: $REF_PROTEOME. Run: bash scripts/download_databases.sh $CONFIG"
    [[ "$REF_NAME_MODE" =~ ^(annotate|fill|prefer)$ ]] \
        || die "REF_NAME_MODE must be annotate, fill or prefer"
fi

mkdir -p "$OUTDIR"/{01_proteins,02_swissprot,03_interproscan,04_eggnog,05_merged,06_final,tmp}
cd "$OUTDIR"

PROT="$OUTDIR/01_proteins/${SPECIES}.longest.prot.fa"
LONGEST_GFF="$OUTDIR/01_proteins/${SPECIES}.longest.gff3"
SPROT_FA="$DB_DIR/uniprot_sprot.fasta"
SPROT_DMND="$DB_DIR/uniprot_sprot.dmnd"
BLAST_OUT="$OUTDIR/02_swissprot/${SPECIES}.sprot.outfmt6.tsv"
IPR_PREFIX="$OUTDIR/03_interproscan/${SPECIES}.ipr"
EMAP_ANNOT="$OUTDIR/04_eggnog/${SPECIES}.emapper.annotations"
MERGED_DIR="$OUTDIR/05_merged/${SPECIES}_agat"
FINAL_GFF="$OUTDIR/06_final/${SPECIES}.functional.gff3"
SUMMARY="$OUTDIR/06_final/${SPECIES}.functional_annotation_summary.tsv"
VERSIONS="$OUTDIR/06_final/${SPECIES}.software_versions.txt"
REF_OUT="$OUTDIR/05_reference/${REF_SPECIES:-none}"
RBH="$REF_OUT/rbh.tsv"
NAME_CHECK="$OUTDIR/06_final/${SPECIES}.reference_name_check.tsv"

[[ -n "${CONDA_SETUP:-}" ]] && eval "$CONDA_SETUP"
eval "$(conda shell.bash hook)"
# One environment holds AGAT, DIAMOND, eggNOG-mapper, Java and InterProScan
conda activate "$CONDA_ENV"

# =============================================================================
# 1. Longest isoform + protein FASTA
# =============================================================================

if [[ ! -s "$LONGEST_GFF" ]]; then
    log "Step 1a: keeping the longest isoform per gene"
    agat_sp_keep_longest_isoform.pl --gff "$GFF" --output "$LONGEST_GFF"
fi

if [[ ! -s "$PROT" ]]; then
    log "Step 1b: extracting proteins"
    agat_sp_extract_sequences.pl --gff "$LONGEST_GFF" --fasta "$GENOME" --protein \
        --output "$OUTDIR/01_proteins/tmp.prot.fa"
    # Drop stop codons (InterProScan rejects '*') and keep only the ID in headers
    sed '/^>/! s/\*//g' "$OUTDIR/01_proteins/tmp.prot.fa" \
        | awk '/^>/{print $1; next}{print}' > "$PROT"
    rm -f "$OUTDIR/01_proteins/tmp.prot.fa"
fi
NPROT=$(grep -c '^>' "$PROT")
log "Proteins (one per gene): $NPROT"

# =============================================================================
# 2. Swiss-Prot names (DIAMOND)
# =============================================================================
if [[ ! -s "$BLAST_OUT" ]]; then
    log "Step 2: DIAMOND blastp vs Swiss-Prot"
    # Rebuild the index if it was made with an incompatible DIAMOND version
    if ! diamond dbinfo -d "$SPROT_DMND" >/dev/null 2>&1; then
        log "  re-indexing Swiss-Prot for this DIAMOND version"
        diamond makedb --in "$SPROT_FA" -d "${SPROT_DMND%.dmnd}" -p "$THREADS"
    fi
    # Standard 12-column outfmt 6, as required by AGAT
    diamond blastp \
        -q "$PROT" -d "$SPROT_DMND" -o "$BLAST_OUT" \
        --outfmt 6 \
        --evalue "$EVALUE" --query-cover "$QCOV" --subject-cover "$SCOV" \
        --max-target-seqs 1 $DIAMOND_MODE \
        -p "$THREADS" --tmpdir "$OUTDIR/tmp"
fi

# =============================================================================
# 3. InterProScan
# =============================================================================
if [[ ! -s "${IPR_PREFIX}.tsv" ]]; then
    log "Step 3: InterProScan"
    [[ -n "${IPR_SETUP:-}" ]] && eval "$IPR_SETUP"

    IPR_OPTS=()
    [[ "$IPR_OFFLINE" == "yes" ]] && IPR_OPTS+=(-dp)
    [[ -n "${IPR_APPL:-}" ]]      && IPR_OPTS+=(-appl "$IPR_APPL")

    # -pa (MetaCyc/Reactome) is deliberately not used: pathways come from
    # eggNOG/KEGG, which suit vertebrates better
    "$IPR_CMD" \
        -i "$PROT" -b "$IPR_PREFIX" \
        -f TSV,GFF3,XML \
        -iprlookup -goterms \
        -cpu "$THREADS" ${IPR_OPTS[@]+"${IPR_OPTS[@]}"} \
        -T "$OUTDIR/tmp/ipr_tmp"
    rm -rf "$OUTDIR/tmp/ipr_tmp"
fi

# =============================================================================
# 4. eggNOG-mapper
# =============================================================================
if [[ ! -s "$EMAP_ANNOT" ]]; then
    log "Step 4: eggNOG-mapper"
    emapper.py \
        -i "$PROT" --itype proteins \
        -m diamond --sensmode "$EMAPPER_SENSMODE" \
        --data_dir "$EGGNOG_DB" \
        --tax_scope "$EMAPPER_TAX_SCOPE" \
        -o "$SPECIES" --output_dir "$OUTDIR/04_eggnog" \
        --temp_dir "$OUTDIR/tmp" \
        --cpu "$THREADS" --override
fi

# =============================================================================
# 5. Merge Swiss-Prot names + InterProScan into the GFF3 (AGAT)
# =============================================================================
if [[ ! -d "$MERGED_DIR" ]]; then
    log "Step 5: merging annotations into the GFF3"
    agat_sp_manage_functional_annotation.pl \
        --gff "$LONGEST_GFF" \
        --blast "$BLAST_OUT" --db "$SPROT_FA" \
        --interpro "${IPR_PREFIX}.tsv" \
        --output "$MERGED_DIR"
fi

AGAT_GFF=$(find "$MERGED_DIR" -maxdepth 1 -name '*.gff*' | head -1)
[[ -s "$AGAT_GFF" ]] || die "no GFF3 found in $MERGED_DIR"

# =============================================================================
# 5b. Optional: reciprocal best hits against a reference species
# =============================================================================
if [[ -n "$REF_SPECIES" && ! -s "$RBH" ]]; then
    log "Step 5b: reciprocal best hits against $REF_SPECIES"
    bash "$PIPE_DIR/scripts/reference_rbh.sh" \
        "$PROT" "$REF_PROTEOME" "$REF_OUT" "$THREADS" "$EVALUE" "$QCOV" "$SCOV"
fi

# =============================================================================
# 6. Post-processing
# =============================================================================
# Rerun when the final GFF3 is missing, or when a reference species was added
if [[ ! -s "$FINAL_GFF" || ( -n "$REF_SPECIES" && ! -s "$NAME_CHECK" ) ]]; then
    log "Step 6: post-processing"
    cp "$AGAT_GFF" "$OUTDIR/tmp/step6.gff3"

    # 6a. Gene-symbol case (only the Name attribute is changed)
    if [[ "$NAME_CASE" == "upper" ]]; then
        sed -E 's/(;Name=)([^;]+)/\1\U\2/' \
            "$OUTDIR/tmp/step6.gff3" \
            | sed -E 's/(Name=C[0-9XYZW]+)ORF([0-9])/\1orf\2/' \
            > "$OUTDIR/tmp/step6a.gff3"
        mv "$OUTDIR/tmp/step6a.gff3" "$OUTDIR/tmp/step6.gff3"
    fi

    # 6b. Reference-species orthologs: compare names, then annotate/fill/prefer
    if [[ -n "$REF_SPECIES" ]]; then
        # Swiss-Prot name vs reference symbol for every reciprocal best hit
        awk -F'\t' -v OFS='\t' '
            FNR == 1 { f++ }
            f == 1 { if (FNR > 1) { rsym[$1] = toupper($5); rid[$1] = $2 }; next }
            $3 == "mRNA" {
                id = $9; sub(/^ID=/, "", id); sub(/;.*/, "", id)
                if (!(id in rid)) next
                nm = ""; if (match($9, /(^|;)Name=[^;]*/)) { nm = substr($9, RSTART, RLENGTH); sub(/^;?Name=/, "", nm) }
                base = toupper(nm); sub(/_[0-9]+$/, "", base)
                cls = (rsym[id] == "-") ? "ref_no_symbol" : (nm == "") ? "no_sprot_name" : (base == rsym[id]) ? "same" : "different"
                print id, rid[id], (nm == "" ? "-" : nm), rsym[id], cls
            }' "$RBH" "$OUTDIR/tmp/step6.gff3" \
            | sort | { echo -e "mRNA_id\tref_gene_id\tsprot_name\tref_symbol\tcomparison"; cat; } > "$NAME_CHECK"

        awk -F'\t' -v OFS='\t' -v mode="$REF_NAME_MODE" -v upper="$NAME_CASE" -v ref="$REF_SPECIES" \
            -f "$PIPE_DIR/scripts/apply_reference_names.awk" \
            "$RBH" "$OUTDIR/tmp/step6.gff3" "$OUTDIR/tmp/step6.gff3" \
            > "$OUTDIR/tmp/step6r.gff3"
        mv "$OUTDIR/tmp/step6r.gff3" "$OUTDIR/tmp/step6.gff3"
    fi

    # 6c. eggNOG Preferred_name for genes still without a name
    if [[ "$EGGNOG_NAME_FALLBACK" == "yes" ]]; then
        awk -F'\t' -v OFS='\t' -v upper="$NAME_CASE" \
            -f "$PIPE_DIR/scripts/add_eggnog_names.awk" \
            "$EMAP_ANNOT" "$OUTDIR/tmp/step6.gff3" "$OUTDIR/tmp/step6.gff3" \
            > "$OUTDIR/tmp/step6b.gff3"
        mv "$OUTDIR/tmp/step6b.gff3" "$OUTDIR/tmp/step6.gff3"
    fi

    # 6d. Strip GO source tags, e.g. GO:0003677(InterPro) -> GO:0003677
    if [[ "$CLEAN_GO_TAGS" == "yes" ]]; then
        sed -E 's/(GO:[0-9]{7})\([A-Za-z]+\)/\1/g' "$OUTDIR/tmp/step6.gff3" \
            > "$OUTDIR/tmp/step6c.gff3"
        mv "$OUTDIR/tmp/step6c.gff3" "$OUTDIR/tmp/step6.gff3"
    fi

    mv "$OUTDIR/tmp/step6.gff3" "$FINAL_GFF"
fi

# =============================================================================
# 7. Summary statistics and software versions
# =============================================================================
log "Step 7: summary"
S="$OUTDIR/tmp/sets"; mkdir -p "$S"
IPR="${IPR_PREFIX}.tsv"

grep '^>' "$PROT" | sed 's/^>//' | sort -u                                 > "$S/all"
cut -f1 "$BLAST_OUT" | sort -u                                             > "$S/swissprot_hit"
awk -F'\t' '$3=="mRNA" && $9 ~ /(^|;)Name=/' "$FINAL_GFF" \
    | grep -o 'ID=[^;]*' | sed 's/ID=//' | sort -u                         > "$S/named_final"
cut -f1 "$IPR" | sort -u                                                   > "$S/interproscan_any"
awk -F'\t' '$12!="" && $12!="-"' "$IPR" | cut -f1 | sort -u                > "$S/interpro_entry"
awk -F'\t' '$4=="Pfam"' "$IPR" | cut -f1 | sort -u                         > "$S/pfam"
awk -F'\t' '$14!="" && $14!="-"' "$IPR" | cut -f1 | sort -u                > "$S/go_interproscan"
grep -v '^#' "$EMAP_ANNOT" | cut -f1 | sort -u                             > "$S/eggnog_hit"
grep -v '^#' "$EMAP_ANNOT" | awk -F'\t' '$10!="-"' | cut -f1 | sort -u     > "$S/go_eggnog"
grep -v '^#' "$EMAP_ANNOT" | awk -F'\t' '$12!="-"' | cut -f1 | sort -u     > "$S/kegg_ko"
cat "$S/go_interproscan" "$S/go_eggnog" | sort -u                          > "$S/go_any"
cat "$S/swissprot_hit" "$S/interpro_entry" "$S/eggnog_hit" | sort -u       > "$S/any_annotation"
CATS="swissprot_hit named_final interproscan_any interpro_entry pfam go_interproscan eggnog_hit go_eggnog kegg_ko go_any any_annotation"
if [[ -n "$REF_SPECIES" ]]; then
    awk -F'\t' 'NR > 1' "$RBH" | cut -f1 | sort -u                          > "$S/ref_ortholog"
    awk -F'\t' 'NR > 1 && $5 == "same"'      "$NAME_CHECK" | cut -f1          > "$S/ref_name_same"
    awk -F'\t' 'NR > 1 && $5 == "different"' "$NAME_CHECK" | cut -f1          > "$S/ref_name_different"
    CATS="$CATS ref_ortholog ref_name_same ref_name_different"
fi

TOTAL=$(wc -l < "$S/all")
{
    printf "category\tn_genes\tpercent\n"
    for f in $CATS; do
        n=$(wc -l < "$S/$f")
        awk -v f="$f" -v n="$n" -v t="$TOTAL" \
            'BEGIN{printf "%s\t%d\t%.1f\n", f, n, 100*n/t}'
    done
    printf "total_genes\t%d\t100.0\n" "$TOTAL"
} > "$SUMMARY"
column -t "$SUMMARY" 2>/dev/null || cat "$SUMMARY"

# Species of the Swiss-Prot entries used for naming (top 15)
grep -P '\tmRNA\t' "$FINAL_GFF" | grep -o 'os=[^;]*' | sed 's/os=//' \
    | sort | uniq -c | sort -rn | head -15 \
    > "$OUTDIR/06_final/${SPECIES}.swissprot_species.txt" || true

# Software versions, for the methods section
{
    echo "Pipeline run: $(date '+%F %T')"
    echo "DIAMOND:       $(diamond version 2>/dev/null | head -1)"
    echo "AGAT:          $(agat_sp_keep_longest_isoform.pl --help 2>/dev/null | grep -o 'Version: v[0-9.]*' | head -1)"
    echo "eggNOG-mapper: $(emapper.py --version --data_dir "$EGGNOG_DB" 2>/dev/null | grep -o 'emapper-[0-9.]*' | head -1)"
    ( [[ -n "${IPR_SETUP:-}" ]] && eval "$IPR_SETUP"
      echo "InterProScan:  $("$IPR_CMD" -version 2>/dev/null | grep -o 'InterProScan version [0-9.-]*' | head -1)" )
    [[ -s "$DB_DIR/reldate.txt" ]] && echo "Swiss-Prot:    $(head -1 "$DB_DIR/reldate.txt")"
    echo "eggNOG DB:     emapperdb-5.0.2"
    [[ -n "$REF_SPECIES" && -s "$REF_DIR/${REF_SPECIES}.source.txt" ]] \
        && echo "Reference:     $(cat "$REF_DIR/${REF_SPECIES}.source.txt")"
} > "$VERSIONS" 2>/dev/null || true

log "Done."
log "Final GFF3:   $FINAL_GFF"
log "eggNOG table: $EMAP_ANNOT"
log "Summary:      $SUMMARY"
