"""Apply optional name policies while preserving GFF attributes and comments."""
import csv
from pathlib import Path
import re
import subprocess
import tempfile


def attributes(value):
    return dict(part.split('=', 1) for part in value.split(';') if '=' in part)


def main(sm):
    text = Path(sm.input.gff).read_text()
    if sm.params.name_case == 'upper':
        lines = []
        for line in text.splitlines(keepends=True):
            fields = line.rstrip('\n').split('\t')
            if not line.startswith('#') and len(fields) == 9:
                fields[8] = ';'.join(
                    'Name=' + item[5:].upper() if item.startswith('Name=') else item
                    for item in fields[8].split(';'))
                line = '\t'.join(fields) + '\n'
            lines.append(line)
        text = ''.join(lines)
    with open(sm.log[0], 'w') as log, tempfile.TemporaryDirectory() as tmp:
        current = Path(tmp) / 'current.gff3'
        current.write_text(text)
        with open(sm.output.comparison, 'w') as out:
            writer = csv.writer(out, delimiter='\t', lineterminator='\n')
            writer.writerow(['mRNA_id', 'ref_gene_id', 'sprot_name', 'ref_symbol', 'comparison'])
            if sm.input.rbh:
                with open(sm.input.rbh[0]) as handle:
                    hits = {row['query_id']: row for row in csv.DictReader(handle, delimiter='\t')}
                comparisons = []
                for line in text.splitlines():
                    fields = line.split('\t')
                    if line.startswith('#') or len(fields) != 9 or fields[2] != 'mRNA':
                        continue
                    attr = attributes(fields[8])
                    hit = hits.get(attr.get('ID'))
                    if hit is None:
                        continue
                    name, symbol = attr.get('Name', ''), hit['ref_symbol']
                    base = re.sub(r'_\d+$', '', name.upper())
                    status = ('ref_no_symbol' if symbol == '-' else 'no_sprot_name' if not name
                              else 'same' if base == symbol.upper() else 'different')
                    comparisons.append([attr['ID'], hit['ref_gene_id'], name or '-', symbol, status])
                writer.writerows(sorted(comparisons))
                result = subprocess.run([
                    'awk', '-F\t', '-v', 'OFS=\t', '-v', f'mode={sm.params.reference_mode}',
                    '-v', f'upper={sm.params.name_case}', '-v', f'ref={sm.params.reference_species}',
                    '-f', str(sm.input.reference_awk), str(sm.input.rbh[0]), str(current), str(current)
                ], check=True, stdout=subprocess.PIPE, stderr=log, text=True)
                current.write_text(result.stdout)
        if sm.params.eggnog_names:
            result = subprocess.run([
                'awk', '-F\t', '-v', 'OFS=\t', '-v', f'upper={sm.params.name_case}',
                '-f', str(sm.input.eggnog_awk), str(sm.input.eggnog), str(current), str(current)
            ], check=True, stdout=subprocess.PIPE, stderr=log, text=True)
            current.write_text(result.stdout)
        text = current.read_text()
        if sm.params.clean_go_tags:
            text = re.sub(r'(GO:\d{7})\([A-Za-z]+\)', r'\1', text)
        Path(sm.output.gff).write_text(text)


main(snakemake)
