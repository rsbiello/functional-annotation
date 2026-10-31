#!/bin/bash
# =============================================================================
# Reciprocal best hits (RBH) between the annotated proteins and the proteome
# of a closely related reference species (Ensembl *.pep.all.fa.gz)
#
# Usage:
#   bash reference_rbh.sh <query_proteins.fa> <reference.pep.all.fa.gz> <outdir> \
#                         <threads> <evalue> <qcov> <scov>
#
# Output: <outdir>/rbh.tsv with columns
#   query_id  ref_gene_id  pident  evalue  ref_symbol  ref_description
# =============================================================================
set -euo pipefail

PROT=$1; REF=$2; OUT=$3; THREADS=$4; EVALUE=$5; QCOV=$6; SCOV=$7
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
mkdir -p "$OUT"

# 1. Longest protein per reference gene
if [[ ! -s "$OUT/ref.fa" ]]; then
    { [[ "$REF" == *.gz ]] && zcat "$REF" || cat "$REF"; } \
        | awk -v tsv="$OUT/ref_genes.tsv" -f "$HERE/prepare_reference.awk" > "$OUT/ref.fa"
fi
echo "Reference genes (longest protein each): $(grep -c '^>' "$OUT/ref.fa")"

# 2. DIAMOND databases
[[ -s "$OUT/ref.dmnd"   ]] || diamond makedb --in "$OUT/ref.fa" -d "$OUT/ref"   -p "$THREADS" --quiet
[[ -s "$OUT/query.dmnd" ]] || diamond makedb --in "$PROT"       -d "$OUT/query" -p "$THREADS" --quiet

COMMON=(--outfmt 6 qseqid sseqid pident evalue bitscore
        --evalue "$EVALUE" --query-cover "$QCOV" --subject-cover "$SCOV"
        --max-target-seqs 1 --sensitive -p "$THREADS" --quiet)

# 3. Forward (query -> reference) and reverse (reference -> query) searches
[[ -s "$OUT/fwd.tsv" ]] || diamond blastp -q "$PROT"       -d "$OUT/ref"   -o "$OUT/fwd.tsv" "${COMMON[@]}"
[[ -s "$OUT/rev.tsv" ]] || diamond blastp -q "$OUT/ref.fa" -d "$OUT/query" -o "$OUT/rev.tsv" "${COMMON[@]}"

# 4. Keep pairs that are each other's best hit; add symbol and description
awk -F'\t' -v OFS='\t' '
    FILENAME == ARGV[1] { if (!($1 in fwd)) { fwd[$1] = $2; pid[$1] = $3; ev[$1] = $4 }; next }
    FILENAME == ARGV[2] { if (!($1 in rev)) rev[$1] = $2; next }
    FILENAME == ARGV[3] { sym[$1] = $2; des[$1] = $3; next }
    END {
        print "query_id", "ref_gene_id", "pident", "evalue", "ref_symbol", "ref_description"
        for (q in fwd) {
            r = fwd[q]
            if (rev[r] == q) print q, r, pid[q], ev[q], (r in sym ? sym[r] : "-"), (r in des ? des[r] : "-")
        }
    }' "$OUT/fwd.tsv" "$OUT/rev.tsv" "$OUT/ref_genes.tsv" \
    | { read -r header; echo "$header"; sort; } > "$OUT/rbh.tsv"

echo "Reciprocal best hits: $(( $(wc -l < "$OUT/rbh.tsv") - 1 ))"
