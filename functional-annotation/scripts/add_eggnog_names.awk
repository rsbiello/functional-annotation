#!/usr/bin/awk -f
# Add eggNOG Preferred_name to mRNAs (and their genes) that have no Name.
# Usage:
#   awk -F'\t' -v OFS='\t' -f add_eggnog_names.awk \
#       emapper.annotations in.gff3 in.gff3 > out.gff3
# The GFF3 is read twice: pass 1 finds unnamed mRNAs, pass 2 writes the file.

function attr(s, key,    re, m) {
    re = "(^|;)" key "=[^;]*"
    if (match(s, re)) {
        m = substr(s, RSTART, RLENGTH)
        sub("^;?" key "=", "", m)
        return m
    }
    return ""
}

FNR == 1 { f++ }

# ---- file 1: eggNOG annotations (col 1 = query, col 9 = Preferred_name)
f == 1 {
    if ($0 !~ /^#/ && $9 != "" && $9 != "-") egg[$1] = (upper == "upper" ? toupper($9) : $9)
    next
}

# ---- file 2: pass 1 over the GFF3
f == 2 {
    if ($0 ~ /^#/ || NF < 9) next
    if ($3 == "mRNA") {
        id = attr($9, "ID"); par = attr($9, "Parent")
        if (attr($9, "Name") == "" && (id in egg)) {
            newname[id] = egg[id]
            if (!(par in genename)) genename[par] = egg[id]
        }
    }
    next
}

# ---- file 3: pass 2, write the output
{
    if ($0 !~ /^#/ && NF >= 9) {
        id = attr($9, "ID")
        if ($3 == "mRNA" && (id in newname)) {
            $9 = $9 ";Name=" newname[id] ";name_source=eggNOG"
            n++
        } else if ($3 == "gene" && (id in genename) && attr($9, "Name") == "") {
            $9 = $9 ";Name=" genename[id]
        }
    }
    print
}

END { printf "eggNOG fallback names added to %d mRNAs\n", n > "/dev/stderr" }
