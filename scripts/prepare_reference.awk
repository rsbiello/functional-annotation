#!/usr/bin/awk -f
# Reduce an Ensembl *.pep.all.fa proteome to the longest protein-coding
# protein per gene.
# Usage:
#   zcat Species.pep.all.fa.gz | awk -v tsv=ref_genes.tsv -f prepare_reference.awk > ref.fa
# Output FASTA headers are Ensembl gene IDs; ref_genes.tsv holds
#   gene_id <tab> gene_symbol <tab> description

function flush() {
    if (gene != "" && (bt == "" || bt == "protein_coding")) {
        gsub(/\*/, "", seq)
        if (length(seq) > len[gene]) {
            len[gene] = length(seq); best[gene] = seq
            sym[gene] = symbol;      des[gene] = desc
        }
    }
    gene = ""; seq = ""
}

/^>/ {
    flush()
    symbol = ""; desc = ""; bt = ""
    for (i = 2; i <= NF; i++) {
        if ($i ~ /^gene:/)                    gene   = substr($i, 6)
        else if ($i ~ /^gene_symbol:/)        symbol = substr($i, 13)
        else if ($i ~ /^transcript_biotype:/) bt     = substr($i, 20)
    }
    if (gene == "") { gene=$1; sub(/^>/,"",gene) }
    else sub(/\.[0-9]+$/, "", gene)               # drop the version suffix
    if (match($0, /description:.*/)) {
        desc = substr($0, RSTART + 12)
        sub(/ *\[Source:.*$/, "", desc)
    }
    next
}

{ seq = seq $0 }

END {
    flush()
    for (g in best) {
        print ">" g
        s = best[g]
        for (i = 1; i <= length(s); i += 60) print substr(s, i, 60)
        print g "\t" (sym[g] == "" ? "-" : sym[g]) "\t" (des[g] == "" ? "-" : des[g]) > tsv
    }
}
