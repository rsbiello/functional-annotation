#!/usr/bin/awk -f
# Add reference-species orthologs (reciprocal best hits) to a GFF3 and,
# depending on the mode, use their gene symbols as names.
# Usage:
#   awk -F'\t' -v OFS='\t' -v mode=annotate|fill|prefer -v upper=upper|keep \
#       -v ref=gallus_gallus -f apply_reference_names.awk rbh.tsv in.gff3 in.gff3 > out.gff3
#
#   annotate : add ref_species, ref_ortholog, ref_symbol, ref_product; keep Name
#   fill     : as annotate, and use ref_symbol as Name where there is no Name
#   prefer   : as annotate, and use ref_symbol as Name whenever it exists
#              (the previous name is kept as sprot_name)

function attr(s, key,    re, m) {
    re = "(^|;)" key "=[^;]*"
    if (match(s, re)) {
        m = substr(s, RSTART, RLENGTH); sub("^;?" key "=", "", m); return m
    }
    return ""
}
function setattr(s, key, val,    re) {
    re = "(^|;)" key "=[^;]*"
    if (match(s, re)) {
        return substr(s, 1, RSTART - 1) (RSTART > 1 ? ";" : "") key "=" val substr(s, RSTART + RLENGTH)
    }
    return s ";" key "=" val
}
function esc(v) {                      # GFF3 reserved characters
    gsub(/%/, "%25", v); gsub(/;/, "%3B", v); gsub(/=/, "%3D", v)
    gsub(/&/, "%26", v); gsub(/,/, "%2C", v); gsub(/\t/, " ", v)
    return v
}

FNR == 1 { if (FILENAME == ARGV[1]) f=1; else if (f < 2) f=2; else f=3 }

# ---- file 1: rbh.tsv
f == 1 {
    if (FNR == 1) next
    rid[$1] = $2
    if ($5 != "-") rsym[$1] = (upper == "upper" ? toupper($5) : $5)
    if ($6 != "-") rdes[$1] = $6
    next
}

# ---- file 2: pass 1, decide new names per mRNA and propagate to genes
f == 2 {
    if ($0 ~ /^#/ || NF < 9 || $3 != "mRNA") next
    id = attr($9, "ID"); par = attr($9, "Parent"); nm = attr($9, "Name")
    if (!(id in rsym)) next
    if (mode == "prefer" || (mode == "fill" && nm == "")) {
        newname[id] = rsym[id]
        if (!(par in gname)) gname[par] = rsym[id]
    }
    next
}

# ---- file 3: pass 2, write the GFF3
{
    if ($0 !~ /^#/ && NF >= 9) {
        id = attr($9, "ID")
        if ($3 == "mRNA" && (id in rid)) {
            $9 = $9 ";ref_species=" ref ";ref_ortholog=" rid[id]
            if (id in rsym) $9 = $9 ";ref_symbol=" esc(rsym[id])
            if (id in rdes) $9 = $9 ";ref_product=" esc(rdes[id])
            if (id in newname) {
                old = attr($9, "Name")
                if (old != "" && old != newname[id]) $9 = $9 ";sprot_name=" old
                $9 = setattr($9, "Name", newname[id])
                $9 = setattr($9, "name_source", "reference")
                n++
            }
        } else if ($3 == "gene" && (id in gname)) {
            $9 = setattr($9, "Name", gname[id])
        }
    }
    print
}

END { printf "Reference names applied to %d mRNAs (mode: %s)\n", n, mode > "/dev/stderr" }
